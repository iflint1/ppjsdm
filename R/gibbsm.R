#' Fit a multivariate Gibbs model to a dataset.
#'
#' @param configuration_list A single configuration or a list of configurations assumed to be drawn from the multivariate Gibbs.
#' @param window Observation window.
#' @param covariates Environmental covariates driving the intensity.
#' @param model String to represent the model we're calibrating. You can check the currently authorised models with a call to `show_short_range_models()`.
#' @param medium_range_model String to represent the model we're calibrating. You can check the currently authorised models with a call to `show_medium_range_models()`.
#' @param short_range Short range interaction radius.
#' @param medium_range Medium range interaction radius.
#' @param long_range Long range interaction radius.
#' @param saturation Saturation parameter of the point process.
#' @param print Print the fitted coefficients?
#' @param use_glmnet Use `glmnet` instead of `glm`?
#' @param use_aic Use AIC instead of BIC for model fitting?
#' @param ndummy Number of dummy points for each type. By default, follows the recommendation of Baddeley et al.
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats AIC BIC binomial coefficients deviance glm.fit lm
#' @importFrom spatstat as.im as.owin
#' @importFrom GA ga
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
                   print = FALSE,
                   use_glmnet = TRUE,
                   use_aic = TRUE,
                   ndummy = 0) {
  parameters <- model_parameters(window = window,
                                 covariates = covariates,
                                 saturation = saturation,
                                 model = model,
                                 medium_range_model = medium_range_model,
                                 short_range = short_range,
                                 medium_range = medium_range,
                                 long_range = long_range)
  short_range <- parameters$short_range
  medium_range <- parameters$medium_range
  long_range <- parameters$long_range
  window <- parameters$window
  covariates <- parameters$covariates
  saturation <- parameters$saturation
  model <- parameters$model
  medium_range_model <- parameters$medium_range_model

  # If we're given a single configuration, convert it to a list.
  if(inherits(configuration_list, "Configuration")) {
    configuration_list <- list(configuration_list)
  }

  # Make sure we're given a list of configurations.
  stopifnot(inherits(configuration_list[[1]], "Configuration"))

  number_configurations <- length(configuration_list)
  number_types <- length(levels(types(configuration_list[[1]])))
  estimate_radii <- is.vector(short_range, mode = "numeric") && length(short_range) == 2
  if(estimate_radii) {
    estimate_alpha <- !all(short_range[1] == short_range[2])
    estimate_gamma <- !all(medium_range[1] == medium_range[2] & long_range[1] == long_range[2] & medium_range[1] == long_range[1])

    lower <- c(rep(short_range[1], number_types + 1), rep(medium_range[1], number_types + 1), rep(long_range[1], number_types + 1))
    upper <- c(rep(short_range[2], number_types + 1), rep(medium_range[2], number_types + 1), rep(long_range[2], number_types + 1))
    initial <- (lower + upper) / 2
    get_fit <- function(v) {
      sh <- diag(v[1:number_types], number_types)
      lower <- lower.tri(sh, diag = FALSE)
      upper <- upper.tri(sh, diag = FALSE)
      sh[lower] <- v[number_types + 1]
      sh[upper] <- t(sh)[upper]

      me <- sh + diag(v[(2 + number_types):(1 + 2 * number_types)], number_types)
      me[lower] <- sh[lower] + v[2 * (1 + number_types)]
      me[upper] <- t(me)[upper]

      lo <- me + diag(v[(3 + 2 * number_types):(2 + 3 * number_types)], number_types)
      lo[lower] <- me[lower] + v[3 * (1 + number_types)]
      lo[upper] <- t(lo)[upper]

      # The fitting procedure samples additional points, let us choose their marks in the same range as current ones.
      mark_range <- c(min(get_marks(configuration_list[[1]])), max(get_marks(configuration_list[[1]])))

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
                                              TRUE,
                                              ndummy = ndummy,
                                              estimate_alpha = estimate_alpha,
                                              estimate_gamma = estimate_gamma)

      fit <- fit_gibbs(gibbsm_data_list, use_glmnet = FALSE, use_aic = use_aic, estimate_alpha = estimate_alpha, estimate_gamma = estimate_gamma)
      list(fit = fit, sh = sh, me = me, lo = lo)
    }

    to_optimise <- function(v) {
      out <- tryCatch(
      {
        fit <- get_fit(v)$fit
        if(use_aic) {
          average <- fit$aic
        } else {
          average <- fit$bic
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
    mark_range <- c(min(get_marks(configuration_list[[1]])), max(get_marks(configuration_list[[1]])))
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
                                            FALSE,
                                            ndummy = ndummy,
                                            estimate_alpha = estimate_alpha,
                                            estimate_gamma = estimate_gamma)

    fitted <- fit_gibbs(gibbsm_data_list, use_glmnet = use_glmnet, use_aic = use_aic, estimate_alpha = estimate_alpha, estimate_gamma = estimate_gamma)
  } else {
    short_range <- as.matrix(short_range)
    medium_range <- as.matrix(medium_range)
    long_range <- as.matrix(long_range)

    estimate_alpha <- !all(short_range == 0)
    estimate_gamma <- !all(medium_range == long_range)

    # The fitting procedure samples additional points, let us choose their marks in the same range as current ones.
    mark_range <- c(min(get_marks(configuration_list[[1]])), max(get_marks(configuration_list[[1]])))
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
                                            FALSE,
                                            ndummy = ndummy,
                                            estimate_alpha = estimate_alpha,
                                            estimate_gamma = estimate_gamma)

    fitted <- fit_gibbs(gibbsm_data_list, use_glmnet = use_glmnet, use_aic = use_aic, estimate_alpha = estimate_alpha, estimate_gamma = estimate_gamma)
  }
  fits <-  fitted$fit
  fits_coefficients <- fitted$coefficients
  aic <- fitted$aic
  bic <- fitted$bic

  if(print) {
    print(fits_coefficients)
  }

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
              fit_algorithm = fitted$fit_algorithm)


  if(estimate_radii) {
    ret <- append(ret, list(coefficients = append(fits_coefficients, list(short_range = best_short, medium_range = best_medium, long_range = best_long))))
  } else {
    ret <- append(ret, list(coefficients = append(fits_coefficients, list(short_range = short_range, medium_range = medium_range, long_range = long_range))))
  }
  class(ret) <- "gibbsm"
  ret
}
