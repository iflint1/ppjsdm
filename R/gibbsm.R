fit_gibbs <- function(gibbsm_data,
                      fitting_package,
                      which_lambda,
                      estimate_alpha,
                      estimate_gamma,
                      types_names,
                      covariates_names,
                      use_regularization,
                      nthreads,
                      refit_glmnet,
                      debug,
                      ...) {
  # Extract a few useful quantities used in most regressions
  regressors <- gibbsm_data$regressors
  if(any(is.na(regressors))) {
    max_indices <- 10
    na_row_indices <- which(rowSums(is.na(regressors)) > 0)[seq_len(max_indices)]
    na_row_indices <- na_row_indices[!is.na(na_row_indices)]
    warning(paste0("Some NA values detected in regression matrix; indices: "),
            na_row_indices)
    print(regressors[na_row_indices, ])
  }

  number_types <- length(gibbsm_data$shift)
  nregressors <- ncol(regressors)

  # Penalty factor used in penalized regressions
  pfactor <- rep(1, nregressors)
  pfactor[startsWith(colnames(regressors), "beta0")] <- 0

  if(fitting_package == "glmnet") {
    # Avoid a bug in glmnet: if intercept = FALSE, and there's a column of ones, it gets ignored by glmnet
    # even though its penalty factor is zero.
    if(all(1 == regressors[, startsWith(colnames(regressors), "beta0")])) {
      regressors[1, startsWith(colnames(regressors), "beta0")] <- 1.00001
    }

    # We don't use an offset explicitely because the call to glmnet above returns nonsensical results or hangs.
    # Instead, we'll use a shift for all the beta0 regressors according to -log(rho).
    arguments <- list(x = Matrix(regressors, sparse = TRUE),
                      y = Matrix(gibbsm_data$response, sparse = TRUE),
                      intercept = FALSE,
                      family = "binomial")
    user_supplied <- list(...)

    if(use_regularization) {
      # Add penalty factor
      arguments <- append(arguments, list(penalty.factor = pfactor))

      # Call glmnet
      arguments <- append(arguments, user_supplied)
      if(debug) {
        tm <- Sys.time()
        message("Calling glmnet...")
      }
      fit <- do.call(glmnet, arguments)
      if(debug) {
        message("Finished call to glmnet.")
        print(Sys.time() - tm)
      }
      if(refit_glmnet > 0) {
        lambda_seq <- fit$lambda
        N <- (1 + refit_glmnet) * length(lambda_seq)
        lambda_new <- vector(mode = 'numeric', length = N)
        lambda_new[1] <- lambda_seq[1]
        by <- (lambda_seq[1] - lambda_seq[2]) / lambda_seq[1]
        for(i in 2:N) {
          lambda_new[i] <- (1 - by) * lambda_new[i - 1]
        }

        arguments$lambda <- lambda_new
        arguments$nlambda <- length(lambda_new)

        if(debug) {
          tm <- Sys.time()
          message("Calling glmnet again...")
        }
        fit <- do.call(glmnet, arguments)
        if(debug) {
          message("Finished call to glmnet.")
          print(Sys.time() - tm)
        }
      }

      # Compute AIC or BIC for the vector of lambdas
      # Reference: https://stackoverflow.com/questions/40920051/r-getting-aic-bic-likelihood-from-glmnet
      tLL <- fit$nulldev - deviance(fit)
      k <- fit$df
      n <- fit$nobs
      aic <- -tLL + 2 * k + 2 * k * (k + 1) / (n - k - 1)
      bic <- log(n) * k - tLL

      # Choose lambda that minimizes AIC or BIC.
      if(which_lambda == "AIC") {
        coef <- coefficients(fit)[, which.min(aic)]
        aic <- min(aic)
        bic <- bic[which.min(aic)]
      } else if(which_lambda == "BIC") {
        coef <- coefficients(fit)[, which.min(bic)]
        aic <- aic[which.min(bic)]
        bic <- min(bic)
      } else if(which_lambda == "smallest") {
        coef <- coefficients(fit)
        coef <- coef[, which.min(fit$lambda)]
        aic <- aic[length(aic)]
        bic <- bic[length(bic)]
      } else {
        stop("Unrecognised option for which_lambda, should be one of 'AIC', 'BIC' or 'smallest'.")
      }
    } else {
      # Call glmnet with a decreasing sequence of lambdas
      if(!("lambda" %in% names(user_supplied))) {
        arguments <- append(arguments, list(lambda = rev(0:99)))
      }
      arguments <- append(arguments, user_supplied)
      if(debug) {
        tm <- Sys.time()
        message("Calling glmnet...")
      }
      fit <- do.call(glmnet, arguments)
      if(debug) {
        message("Finished call to glmnet.")
        print(Sys.time() - tm)
      }

      # Compute AIC or BIC for the vector of lambdas
      # Reference: https://stackoverflow.com/questions/40920051/r-getting-aic-bic-likelihood-from-glmnet
      tLL <- fit$nulldev - deviance(fit)
      k <- fit$df
      n <- fit$nobs
      aic <- -tLL + 2 * k + 2 * k * (k + 1) / (n - k - 1)
      bic <- log(n) * k - tLL

      # Choose AIC / BIC corresponding to lambda = 0
      aic <- aic[length(aic)]
      bic <- bic[length(bic)]

      # Extract the corresponding coefficients
      coef <- coefficients(fit, s = 0)
      coef <- setNames(as.vector(coef), nm = rownames(coef))
    }

    # Even with intercept = FALSE, glmnet still returns an intercept 0-valued coefficient
    coef <- coef[-which(names(coef) == "(Intercept)")]

    # We didn't actually use an offset so we have to shift back the intercept coefficients
    shift <- gibbsm_data$shift
    for(i in seq_len(number_types)) {
      coef[match(paste0("beta0_", i), names(coef))] <- coef[match(paste0("beta0_", i), names(coef))] - shift[i]
    }

    # Clean up
    fit_algorithm <- "glmnet"
  } else if(fitting_package == "oem") {

    arguments <- list(x = regressors,
                      y = gibbsm_data$response,
                      intercept = FALSE,
                      family = "binomial",
                      ncores = nthreads)
    user_supplied <- list(...)

    if(!("compute.loss" %in% names(user_supplied))) {
      arguments <- append(arguments, list(compute.loss = TRUE))
    }

    if(use_regularization) {
      # Add penalty factor
      arguments <- append(arguments, list(penalty.factor = pfactor))

      # Add penalty type
      if(!("penalty" %in% names(user_supplied))) {
        arguments <- append(arguments, list(penalty = "lasso"))
      }

      # Call oem
      arguments <- append(arguments, user_supplied)
      if(debug) {
        tm <- Sys.time()
        message("Calling oem...")
      }
      fit <- do.call(oem, arguments)
      if(debug) {
        message("Finished call to oem.")
        print(Sys.time() - tm)
      }

      # Extract coefficients
      coef <- fit$beta$lasso

      # Compute AIC/BIC for vector of lambdas
      # Reference: https://stackoverflow.com/questions/40920051/r-getting-aic-bic-likelihood-from-glmnet
      k <- sapply(seq_len(ncol(coef)), function(col) sum(coef[, col] != 0.))
      n <- fit$nobs
      LL <- 2 * logLik(fit)
      aic <- -LL + 2 * k + 2 * k * (k + 1) / (n - k - 1)
      bic <- -LL + log(n) * k

      # Choose lambda that minimizes AIC or BIC
      if(which_lambda == "AIC") {
        coef <- coef[, which.min(aic)]
        aic <- min(aic)
        bic <- bic[which.min(aic)]
      } else if(which_lambda == "BIC") {
        coef <- coef[, which.min(bic)]
        aic <- aic[which.min(bic)]
        bic <- min(bic)
      } else if(which_lambda == "smallest") {
        coef <- coef[, ncol(coef)]
        aic <- aic[length(aic)]
        bic <- bic[length(bic)]
      } else {
        stop("Unrecognised option for which_lambda, should be one of 'AIC', 'BIC' or 'smallest'.")
      }

      # Put coef in correct format
      coef <- setNames(as.vector(coef), nm = names(coef))
    } else {
      # Use 'ols' penalty for non-regularized regression
      if(!("penalty" %in% names(user_supplied))) {
        arguments <- append(arguments, list(penalty = "ols"))
      }

      # Call oem
      arguments <- append(arguments, user_supplied)
      if(debug) {
        tm <- Sys.time()
        message("Calling oem...")
      }
      fit <- do.call(oem, arguments)
      if(debug) {
        message("Finished call to oem.")
        print(Sys.time() - tm)
      }

      # Extract AIC/BIC using oem functions
      aic <- AIC(fit)
      bic <- BIC(fit)

      # Extract coefficients
      coef <- fit$beta$ols
      coef <- setNames(as.vector(coef), nm = rownames(coef))
    }

    # Remove intercept
    coef <- coef[-which(names(coef) == "(Intercept)")]

    # We didn't actually use an offset so we have to shift back the intercept coefficients
    shift <- gibbsm_data$shift
    for(i in seq_len(number_types)) {
      coef[match(paste0("beta0_", i), names(coef))] <- coef[match(paste0("beta0_", i), names(coef))] - shift[i]
    }

    # Clean up
    fit_algorithm <- "oem"
  } else if(fitting_package == "h2o") {
    # Initialization
    h2o.init(nthreads = nthreads)

    # Put everything into a unique data.frame
    regressors <- as.data.frame(regressors)
    regressors$response <- gibbsm_data$response
    regressors$offset <- gibbsm_data$offset

    arguments <- list(training_frame = as.h2o(regressors),
                      y = 'response',
                      intercept = FALSE,
                      family = "binomial",
                      offset_column = 'offset',
                      standardize = FALSE,
                      link = 'logit')
    user_supplied <- list(...)

    if(use_regularization) {
      stop("Not implemented yet")
    } else {
      # Non-regularized corresponds to lambda = 0
      if(!("lambda" %in% names(user_supplied))) {
        arguments <- append(arguments, list(lambda = 0))
      }

      # Call h2o
      arguments <- append(arguments, user_supplied)
      if(debug) {
        tm <- Sys.time()
        message("Calling h2o")
      }
      fit <- do.call(h2o.glm, arguments)
      if(debug) {
        message("Finished call to h2o")
        print(Sys.time() - tm)
      }

      # Compute AIC/BIC
      warning("AIC/BIC not implemented yet for h2o fitting package")
      aic <- NA
      bic <- NA

      # Extract coefs
      coef <- h2o.coef(fit)
    }

    # Remove intercept from coefficients
    coef <- coef[-which(names(coef) == "Intercept")]

    # Shift intercepts back
    shift <- gibbsm_data$shift
    for(i in seq_len(number_types)) {
      coef[match(paste0("beta0_", i), names(coef))] <- coef[match(paste0("beta0_", i), names(coef))] - shift[i]
    }

    # Clean up
    fit_algorithm <- "h2o"
  } else if(fitting_package == "glm") {
    if(use_regularization) {
      warning("Cannot have regularization with glm, doing a regular glm fit.")
    }

    # Formula for the GLM regression
    fmla <- as.formula(paste0("response ~ 0 + offset(offset) + ", paste0(colnames(regressors), collapse = ' + ')))

    # Put all the required vectors into a data.frame to send to glm
    regressors <- as.data.frame(regressors)
    regressors$offset <- gibbsm_data$offset
    regressors$response <- gibbsm_data$response

    # Call GLM
    if(debug) {
      tm <- Sys.time()
      message("Calling glm...")
    }
    fit <- glm(fmla, family = binomial(), data = regressors, ...)
    if(debug) {
      message("Finished call to glm.")
      print(Sys.time() - tm)
    }

    # Compute AIC/BIC
    aic <- AIC(fit)
    bic <- BIC(fit)

    # Clean up
    coef <- coefficients(fit)
    fit_algorithm <- "glm"
  } else {
    stop("Supplied fitting package not recognized.")
  }

  # Construct coefficients in matrix form
  beta0 <- vector(mode = "numeric", length = number_types)
  alpha <- matrix(0, number_types, number_types)
  gamma <- matrix(0, number_types, number_types)

  # Environmental covariates in vector form
  beta_vector <- coef[!(startsWith(names(coef), "beta0") | startsWith(names(coef), "alpha") | startsWith(names(coef), "gamma") | ("(Intercept)" == names(coef)))]

  # Environmental covariates matrix
  beta <- matrix(NA, nrow = number_types, ncol = length(beta_vector) / number_types)

  # Fill all the coefficients
  for(i in seq_len(number_types)) {
    beta0[i] <- coef[match(paste0("beta0_", i), names(coef))]
    if(length(beta_vector) != 0) {
      beta[i, ] <- beta_vector[endsWith(names(beta_vector), paste0("_", i))]
    }
    if(estimate_alpha[i, i]) {
      alpha[i, i] <- coef[match(paste0("alpha_", i, "_", i), names(coef))]
    }
    if(estimate_gamma[i, i]) {
      gamma[i, i] <- coef[match(paste0("gamma_", i, "_", i), names(coef))]
    }
    if(i < number_types) {
      for(j in (i + 1):number_types) {
        if(estimate_alpha[i, j]) {
          alpha[i, j] <- alpha[j, i] <- coef[match(paste0("alpha_", i, "_", j), names(coef))]
        }
        if(estimate_gamma[i, j]) {
          gamma[i, j] <- gamma[j, i] <- coef[match(paste0("gamma_", i, "_", j), names(coef))]
        }
      }
    }
  }

  # Fill names
  if(!missing(types_names)) {
    names(beta0) <- types_names
    rownames(beta) <- types_names
    rownames(alpha) <- types_names
    colnames(alpha) <- types_names
    rownames(gamma) <- types_names
    colnames(gamma) <- types_names
  }
  if(!missing(covariates_names)) {
    colnames(beta) <- covariates_names
  }

  # Clean up
  list(fit = fit,
       coefficients = list(beta0 = beta0,
                           alpha = alpha,
                           gamma = gamma,
                           beta = beta),
       coefficients_vector = coef,
       aic = aic,
       bic = bic,
       fit_algorithm = fit_algorithm)
}


#' Fit a multivariate Gibbs model to a dataset.
#'
#' @param configuration_list A single configuration or a list of configurations assumed to be drawn from the multivariate Gibbs.
#' @param window Observation window.
#' @param covariates Environmental covariates driving the intensity.
#' @param model String to represent the model we're calibrating. You can check the currently authorised models with a call to show_short_range_models().
#' @param medium_range_model String to represent the model we're calibrating. You can check the currently authorised models with a call to show_medium_range_models().
#' @param short_range Short range interaction radius. Filled with 0.1 by default.
#' @param medium_range Medium range interaction radius. Filled with 0 by default.
#' @param long_range Long range interaction radius. Filled with 0 by default.
#' @param saturation Saturation parameter of the point process.
#' @param debug Print debugging information?
#' @param fitting_package Choices are glmnet, oem or glm.
#' @param which_lambda Which lambda to choose in the regularised fit? Choices are "AIC" (lambda that minimises AIC), "BIC" (lambda that minimises BIC) or "smallest" (smallest lambda).
#' @param max_dummy Maximum number of dummy points for each type. By default, follows the recommendation of Baddeley et al.
#' @param min_dummy Minimum number of dummy points for each type. By default, follows the recommendation of Baddeley et al.
#' @param dummy_factor How many more dummy points there are, compared to the points in the configuration. By default, follows the recommendation of Baddeley et al.
#' @param nthreads Maximum number of threads for parallel computing.
#' @param use_regularization Use the fitting package without regularization.
#' @param return_raw_fitting_data Return the raw fitting data, before calling the GLM fitting package. Mostly used for debugging purposes on large datasets.
#' @param refit_glmnet How many times more lambdas to use on a re-run of glmnet (currently only used with fitting_package = glmnet). For example, refit_glmnet = 0.5 re-runs the fitting procedure with a 50\% longer lambda sequence.
#' @param dummy_distribution How should the dummy distribution be drawn? Choices are "binomial" or "stratified".
#' @param ... Forwarded to the fitting package.
#' @importFrom GA ga
#' @importFrom glmnet glmnet
#' @importFrom h2o as.h2o h2o.coef h2o.glm h2o.init
#' @importFrom Matrix Matrix
#' @importFrom oem oem
#' @importFrom stats AIC as.formula BIC binomial coefficients deviance glm lm logLik setNames
#' @importFrom spatstat.geom as.im as.owin
#' @export
gibbsm <- function(configuration_list,
                   window,
                   covariates,
                   model,
                   medium_range_model,
                   short_range,
                   medium_range,
                   long_range,
                   saturation,
                   debug = FALSE,
                   fitting_package = "glmnet",
                   which_lambda = "AIC",
                   max_dummy = 0,
                   min_dummy = 0,
                   dummy_factor = 0,
                   nthreads = 4,
                   use_regularization = FALSE,
                   return_raw_fitting_data = FALSE,
                   refit_glmnet = 0,
                   dummy_distribution = "binomial",
                   ...) {
  # refit_glmnet only has an effect with glmnet at the moment, do a quick check.
  if(refit_glmnet > 0 & fitting_package != "glmnet") {
    warning("refit_glmnet currently has no effect if the fitting package is not glmnet.")
  }

  # If we're given a single configuration, convert it to a list.
  if(inherits(configuration_list, "Configuration")) {
    configuration_list <- list(configuration_list)
  }

  if(!missing(short_range)) {
    estimate_radii <- is.vector(short_range, mode = "numeric") && length(short_range) == 2
  } else {
    estimate_radii <- FALSE
  }
  parameters <- model_parameters(window = window,
                                 covariates = covariates,
                                 saturation = saturation,
                                 model = model,
                                 medium_range_model = medium_range_model,
                                 short_range = short_range,
                                 medium_range = medium_range,
                                 long_range = long_range,
                                 default_number_types = nlevels(types(configuration_list[[1]])))
  if(!estimate_radii) {
    short_range <- parameters$short_range
    medium_range <- parameters$medium_range
    long_range <- parameters$long_range
  }
  window <- parameters$window
  covariates <- parameters$covariates
  saturation <- parameters$saturation
  model <- parameters$model
  medium_range_model <- parameters$medium_range_model

  # Try to force conversion to Configuration class.
  lapply(configuration_list, function(configuration) as.Configuration(configuration))

  # Set types names and covariates names
  types_names <- levels(types(configuration_list[[1]]))
  covariates_names <- names(covariates)

  # Make sure we're given a list of configurations.
  stopifnot(inherits(configuration_list[[1]], "Configuration"))

  number_configurations <- length(configuration_list)
  number_types <- length(levels(types(configuration_list[[1]])))
  if(estimate_radii) {
    estimate_alpha <- matrix(short_range[1] != short_range[2], nrow = number_types, ncol = number_types)
    estimate_gamma <- matrix(long_range[1] != long_range[2], nrow = number_types, ncol = number_types)

    lower <- c(rep(short_range[1], number_types * (number_types + 1) / 2), rep(medium_range[1], number_types * (number_types + 1) / 2), rep(long_range[1], number_types * (number_types + 1) / 2))
    upper <- c(rep(short_range[2], number_types * (number_types + 1) / 2), rep(medium_range[2], number_types * (number_types + 1) / 2), rep(long_range[2], number_types * (number_types + 1) / 2))
    initial <- (lower + upper) / 2
    get_fit <- function(v) {
      sh <- matrix(NA, number_types, number_types)
      index <- 1
      for(i in seq_len(number_types)) {
        for(j in i:number_types) {
          sh[i, j] <- sh[j, i] <- v[index]
          index <- index + 1
        }
      }

      me <- matrix(NA, number_types, number_types)
      for(i in seq_len(number_types)) {
        for(j in i:number_types) {
          me[i, j] <- me[j, i] <- sh[i, j] + v[index]
          index <- index + 1
        }
      }

      lo <- matrix(NA, number_types, number_types)
      for(i in seq_len(number_types)) {
        for(j in i:number_types) {
          lo[i, j] <- lo[j, i] <- me[i, j] + v[index]
          index <- index + 1
        }
      }

      # The fitting procedure samples additional points, let us choose their marks in the same range as current ones.
      mark_range <- c(min(marks(configuration_list[[1]])), max(marks(configuration_list[[1]])))

      gibbsm_data_list <- prepare_gibbsm_data(configuration_list,
                                              window,
                                              covariates,
                                              model,
                                              medium_range_model,
                                              sh,
                                              me,
                                              lo,
                                              saturation,
                                              mark_range,
                                              max_dummy = max_dummy,
                                              min_dummy = min_dummy,
                                              dummy_factor = dummy_factor,
                                              estimate_alpha = estimate_alpha,
                                              estimate_gamma = estimate_gamma,
                                              nthreads = nthreads,
                                              debug = debug,
                                              dummy_distribution = dummy_distribution)

      fit <- fit_gibbs(gibbsm_data_list,
                       fitting_package = "glm",
                       which_lambda = which_lambda,
                       estimate_alpha = estimate_alpha,
                       estimate_gamma = estimate_gamma,
                       types_names = types_names,
                       covariates_names = covariates_names,
                       use_regularization = use_regularization,
                       nthreads = nthreads,
                       refit_glmnet = refit_glmnet,
                       debug = debug,
                       ...)
      list(fit = fit, sh = sh, me = me, lo = lo)
    }

    to_optimise <- function(v) {
      out <- tryCatch(
      {
        fit <- get_fit(v)$fit
        if(which_lambda == "BIC") {
          average <- fit$bic
        } else {
          average <- fit$aic
        }
        average
      },
      error = function(cond) {
        return(Inf)
      },
      warning = function(cond) {
        return(Inf)
      })
      out
    }
    # MC
    # best_metric <- Inf
    # for(i in seq_len(100)) {
    #   v <- runif(3 * (1 + number_types), lower, upper)
    #   metric <- to_optimise(v)
    #   if(metric < best_metric) {
    #     best_metric <- metric
    #     best_v <- v
    #   }
    # }
    # result <- get_fit(best_v)

    # GA
    GA <- ga(type = "real-valued",
             fitness =  function(v) -to_optimise(v),
             lower = lower,
             upper = upper)
    result <- get_fit(GA@solution)

    # optim
    # opt <- optim(initial, to_optimise, lower = lower, upper = upper, method = "L-BFGS-B")
    # result <- get_fit(opt$par)

    best_short <- result$sh
    best_medium <- result$me
    best_long <- result$lo

    # Refit with best radii and non-approximation
    # The fitting procedure samples additional points, let us choose their marks in the same range as current ones.
    mark_range <- c(min(marks(configuration_list[[1]])), max(marks(configuration_list[[1]])))
    gibbsm_data_list <- prepare_gibbsm_data(configuration_list,
                                            window,
                                            covariates,
                                            model,
                                            medium_range_model,
                                            best_short,
                                            best_medium,
                                            best_long,
                                            saturation,
                                            mark_range,
                                            max_dummy = max_dummy,
                                            min_dummy = min_dummy,
                                            dummy_factor = dummy_factor,
                                            estimate_alpha = estimate_alpha,
                                            estimate_gamma = estimate_gamma,
                                            nthreads = nthreads,
                                            debug = debug,
                                            dummy_distribution = dummy_distribution,
                                            type_names = types_names)

    if(return_raw_fitting_data) {
      return(gibbsm_data_list)
    }

    fitted <- fit_gibbs(gibbsm_data_list,
                        fitting_package = fitting_package,
                        which_lambda = which_lambda,
                        estimate_alpha = estimate_alpha,
                        estimate_gamma = estimate_gamma,
                        types_names = types_names,
                        covariates_names = covariates_names,
                        use_regularization = use_regularization,
                        nthreads = nthreads,
                        refit_glmnet = refit_glmnet,
                        debug = debug,
                        ...)
  } else {
    short_range <- as.matrix(short_range)
    medium_range <- as.matrix(medium_range)
    long_range <- as.matrix(long_range)

    if(nrow(short_range) != number_types |
       ncol(short_range) != number_types |
       nrow(medium_range) != number_types |
       ncol(medium_range) != number_types |
       nrow(long_range) != number_types |
       ncol(long_range) != number_types) {
      stop("The interaction radii have not been provided in the right format.")
    }


    estimate_alpha <- short_range != 0
    estimate_gamma <- medium_range != long_range

    # The fitting procedure samples additional points, let us choose their marks in the same range as current ones.
    mark_range <- c(min(marks(configuration_list[[1]])), max(marks(configuration_list[[1]])))
    gibbsm_data_list <- prepare_gibbsm_data(configuration_list,
                                            window,
                                            covariates,
                                            model,
                                            medium_range_model,
                                            short_range,
                                            medium_range,
                                            long_range,
                                            saturation,
                                            mark_range,
                                            max_dummy = max_dummy,
                                            min_dummy = min_dummy,
                                            dummy_factor = dummy_factor,
                                            estimate_alpha = estimate_alpha,
                                            estimate_gamma = estimate_gamma,
                                            nthreads = nthreads,
                                            debug = debug,
                                            dummy_distribution = dummy_distribution,
                                            type_names = types_names)

    if(return_raw_fitting_data) {
      return(gibbsm_data_list)
    }

    fitted <- fit_gibbs(gibbsm_data_list,
                        fitting_package = fitting_package,
                        which_lambda = which_lambda,
                        estimate_alpha = estimate_alpha,
                        estimate_gamma = estimate_gamma,
                        types_names = types_names,
                        covariates_names = covariates_names,
                        use_regularization = use_regularization,
                        nthreads = nthreads,
                        refit_glmnet = refit_glmnet,
                        debug = debug,
                        ...)
  }
  fits <-  fitted$fit
  fits_coefficients <- fitted$coefficients
  aic <- fitted$aic
  bic <- fitted$bic

  ret <- list(complete = fits,
              configuration_list = configuration_list,
              estimate_alpha = estimate_alpha,
              estimate_gamma = estimate_gamma,
              data_list = gibbsm_data_list,
              parameters = list(model = model,
                                medium_range_model = medium_range_model,
                                covariates = covariates,
                                saturation = saturation),
              coefficients_vector = fitted$coefficients_vector,
              aic = aic,
              bic = bic,
              window = window,
              fit_algorithm = fitted$fit_algorithm,
              used_regularization = use_regularization,
              debug = debug,
              nthreads = nthreads,
              dummy_distribution = dummy_distribution,
              type_names = types_names,
              mark_range = mark_range)

  if(estimate_radii) {
    short_range <- best_short
    medium_range <- best_medium
    long_range <- best_long
  }
  rownames(short_range) <- types_names
  colnames(short_range) <- types_names
  rownames(medium_range) <- types_names
  colnames(medium_range) <- types_names
  rownames(long_range) <- types_names
  colnames(long_range) <- types_names

  ret <- append(ret, list(coefficients = append(fits_coefficients,
                                                list(short_range = short_range,
                                                     medium_range = medium_range,
                                                     long_range = long_range))))
  class(ret) <- "gibbsm"
  ret
}
