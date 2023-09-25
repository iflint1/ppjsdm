#' Summary of a Fitted Gibbs point process.
#' The function can be supplied multiple fits, in which case the aggregate estimator is computed.
#'
#' @param object A fitted model object.
#' @param ... Other fitted model objects.
#' @param list List of fits to consider, in addition to object.
#' @param debug Display debug information?
#' @param time_limit Time limit measured in  `unit` that can be spent running this function.
#' @param unit Unit used to measure the time limit (hours, mins, secs, etc).
#' @param nthreads (optional) number of threads to use.
#' @param npoints Target number of points in the restricted window that the vcov matrix is computed on. Computation is slower for larger values, but the vcov matrix is then better approximated.
#' @param multiple_windows Compute A2 and A3 on a lot of small windows and which are then averaged out, or only on a single restricted window?
#' @param assume_fitted_to_same_data Should the function assume that the model was fitted to the same data?
#' If so and if multiple fits are supplied, the function assumes that all configurations are identical.
#' If not, the function assumes that the different fits were obtained through different realisations of the same point process.
#' @importFrom stats pnorm qnorm
#' @export
summary.gibbsm <- function(object,
                           ...,
                           list,
                           debug = FALSE,
                           time_limit = Inf,
                           unit = "hours",
                           nthreads = NULL,
                           npoints = 2000,
                           multiple_windows = TRUE,
                           assume_fitted_to_same_data = FALSE) {
  # time_limit below is the time_limit to run each of the 4 matrix constructions,
  # so allow time_limit / 4 for each one of them.
  time_limit <- time_limit / 4

  # Allow for either sequence of fits or list of fits, convert both to list
  fits <- if(missing(list)) {
    if(missing(object)) {
      base::list(...)
    } else {
      base::list(object, ...)
    }
  } else {
    if(missing(object)) {
      list
    } else {
      c(list(object), list)
    }
  }

  # Make sure thetas are compatible
  theta1 <- fits[[1]]$coefficients_vector
  for(fit in fits) {
    if(length(theta1) != length(fit$coefficients_vector)) {
      stop("Thetas of the supplied fits do not have the same length.")
    }
    if(!identical(names(theta1), names(fit$coefficients_vector))) {
      stop("Thetas of the supplied fits do not have the same names.")
    }
  }

  # Make sure estimate_alpha and estimate_gamma are all identical
  estimate_alpha <- fits[[1]]$estimate_alpha
  estimate_gamma <- fits[[1]]$estimate_gamma
  for(fit in fits) {
    for(i in seq_len(length(estimate_alpha))) {
      if(!all(estimate_alpha[[i]] == fit$estimate_alpha[[i]])) {
        stop("estimate_alpha is not the same for all fits.")
      }
    }
    if(!all(estimate_gamma == fit$estimate_gamma)) {
      stop("estimate_gamma is not the same for all fits.")
    }
  }

  # Make sure compatible options were used
  for(fit in fits) {
    if(fit$used_regularization) {
      warning("Computing the Variance-Covariance matrix of a regularised fit.")
    }
    if(length(fit$configuration_list) != 1) {
      stop("Cannot compute VCOV matrix for a fit obtained on a list of configurations.")
    }
  }

  # Compute the regression coefficient, averaged out over the fits
  average_theta <- sapply(names(theta1), function(nm) {
    mean(sapply(fits, function(fit) fit$coefficients_vector[nm]), na.rm = TRUE)
  })

  # Compute fixed parameters
  number_types <- length(fits[[1]]$coefficients$beta0)
  types_names <- names(fits[[1]]$coefficients$beta0)
  covariates_names <- names(fits[[1]]$covariates)

  # Make sure all fits have the same configuration
  if(!assume_fitted_to_same_data) {
    configuration <- fits[[1]]$configuration_list[[1]]
    for(fit in fits) {
      if(!identical(configuration, fit$configuration_list[[1]])) {
        if(length(configuration$x) != length(fit$configuration_list[[1]]$x)) {
          # Found two configurations with different number of points, assuming that the data was fitted to different data.
          assume_fitted_to_same_data <- FALSE
        } else if(!identical(types(configuration), types(fit$configuration_list[[1]])) |
                  !identical(marks(configuration), marks(fit$configuration_list[[1]]))) {
          # Found two configurations with same number of points, but different types/marks. Assuming that the data was fitted to different data.
          assume_fitted_to_same_data <- FALSE
        } else {
          error <- sqrt(0.5 * (mean((configuration$x - fit$configuration_list[[1]]$x)^2) +
                                 mean((configuration$y - fit$configuration_list[[1]]$y)^2)))
          if(error < 1e4) {
            # Found two configurations with same number of points and same types/marks, and difference between locations.
            # The locations were probably jittered around.
            assume_fitted_to_same_data <- TRUE
          } else {
            assume_fitted_to_same_data <- FALSE
          }
        }
      }
    }
  }

  # Start constructing the return object
  y <- base::list()
  class(y) <- "summary_gibbsm"

  S <- compute_S(list = fits, debug = debug, time_limit = time_limit, unit = unit)
  A1 <- compute_A1(list = fits, nthreads = nthreads, debug = debug,
                   time_limit = time_limit, unit = unit)
  A2A3 <- compute_A2_plus_A3(list = fits, nthreads = nthreads, debug = debug,
                                   npoints = npoints, multiple_windows = multiple_windows,
                                   time_limit = time_limit, unit = unit)
  G2 <- compute_G2(list = fits, nthreads = nthreads, debug = debug,
                   time_limit = time_limit, unit = unit)

  S_inv <- solve(S)
  SA1S <- S_inv %*% A1 %*% S_inv
  SA2A3S <- S_inv %*% A2A3 %*% S_inv
  SG2S <- S_inv %*% G2 %*% S_inv

  if(!assume_fitted_to_same_data) {
    SA1S <- SA1S / length(fits)
    SA2A3S <- SA2A3S / length(fits)
  }

  se_numerical <- sqrt(diag(SG2S))
  se <- sqrt(diag(SA1S) + diag(SA2A3S) + diag(SG2S))

  se_numerical_proportion <- se_numerical / se
  coefficients <- average_theta
  one_ninetysix <- qnorm(0.975)
  lo <- coefficients - one_ninetysix * se
  hi <- coefficients + one_ninetysix * se
  lo_numerical <- coefficients - one_ninetysix * se_numerical
  hi_numerical <- coefficients + one_ninetysix * se_numerical
  zval <- coefficients / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  psig <- cut(pval,
              breaks = c(0, 0.001, 0.01, 0.05, 1),
              labels = c("***", "**", "*", "  "),
              include.lowest = TRUE)
  y$G1 <- SA1S + SA2A3S
  y$G2 <- SG2S
  y$coefficients <- data.frame(coefficients = coefficients,
                               se = se,
                               CI95_lo = lo,
                               CI95_hi = hi,
                               Ztest = psig,
                               Pval = pval,
                               Zval = zval,
                               se_numerical_proportion = se_numerical_proportion)

  y$se <- ppjsdm::format_coefficient_vector(coefficient_vector = se,
                                            number_types = number_types,
                                            types_names = types_names,
                                            covariates_names = covariates_names,
                                            estimate_alpha = estimate_alpha,
                                            estimate_gamma = estimate_gamma)
  y$se_numerical <- ppjsdm::format_coefficient_vector(coefficient_vector = se_numerical,
                                                      number_types = number_types,
                                                      types_names = types_names,
                                                      covariates_names = covariates_names,
                                                      estimate_alpha = estimate_alpha,
                                                      estimate_gamma = estimate_gamma)
  y$lo <- ppjsdm::format_coefficient_vector(coefficient_vector = lo,
                                            number_types = number_types,
                                            types_names = types_names,
                                            covariates_names = covariates_names,
                                            estimate_alpha = estimate_alpha,
                                            estimate_gamma = estimate_gamma)
  y$hi <- ppjsdm::format_coefficient_vector(coefficient_vector = hi,
                                            number_types = number_types,
                                            types_names = types_names,
                                            covariates_names = covariates_names,
                                            estimate_alpha = estimate_alpha,
                                            estimate_gamma = estimate_gamma)
  y$lo_numerical <- ppjsdm::format_coefficient_vector(coefficient_vector = lo_numerical,
                                                      number_types = number_types,
                                                      types_names = types_names,
                                                      covariates_names = covariates_names,
                                                      estimate_alpha = estimate_alpha,
                                                      estimate_gamma = estimate_gamma)
  y$hi_numerical <- ppjsdm::format_coefficient_vector(coefficient_vector = hi_numerical,
                                                      number_types = number_types,
                                                      types_names = types_names,
                                                      covariates_names = covariates_names,
                                                      estimate_alpha = estimate_alpha,
                                                      estimate_gamma = estimate_gamma)

  y

  # vc <- compute_vcov(configuration = object$configuration_list[[1]],
  #                    dummy = object$data_list$dummy,
  #                    window = object$window,
  #                    covariates = object$parameters$covariates,
  #                    model = object$parameters$model,
  #                    medium_range_model = object$parameters$medium_range_model,
  #                    short_range = object$coefficients$short_range,
  #                    medium_range = object$coefficients$medium_range,
  #                    long_range = object$coefficients$long_range,
  #                    saturation = object$parameters$saturation,
  #                    rho = exp(-object$data_list$shift),
  #                    theta = object$coefficients_vector,
  #                    regressors = as.matrix(object$data_list$regressors),
  #                    data_list = object$data_list,
  #                    estimate_alpha = object$estimate_alpha,
  #                    estimate_gamma = object$estimate_gamma,
  #                    debug = object$debug,
  #                    nthreads = object$nthreads,
  #                    npoints = npoints,
  #                    multiple_windows = multiple_windows,
  #                    dummy_distribution = object$dummy_distribution,
  #                    mark_range = object$mark_range)
  #
  # se_numerical <- sqrt(diag(vc$G2))
  # se <- sqrt(diag(vc$G1 + vc$G2))
  #
  # se_numerical_proportion <- se_numerical / se
  # coefficients <- object$coefficients_vector
  # one_ninetysix <- qnorm(0.975)
  # lo <- coefficients - one_ninetysix * se
  # hi <- coefficients + one_ninetysix * se
  # lo_numerical <- coefficients - one_ninetysix * se_numerical
  # hi_numerical <- coefficients + one_ninetysix * se_numerical
  # zval <- coefficients / se
  # pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  # psig <- cut(pval,
  #             breaks = c(0,0.001, 0.01, 0.05, 1),
  #             labels = c("***", "**", "*", "  "),
  #             include.lowest = TRUE)
  # y$G1 <- vc$G1
  # y$G2 <- vc$G2
  # y$coefficients <- data.frame(coefficients = coefficients,
  #                              se = se,
  #                              CI95_lo = lo,
  #                              CI95_hi = hi,
  #                              Ztest = psig,
  #                              Pval = pval,
  #                              Zval = zval,
  #                              se_numerical_proportion = se_numerical_proportion)
  #
  # y$se <- format_coefficient_vector(coefficient_vector = se,
  #                                   number_types = nlevels(types(object$configuration_list[[1]])),
  #                                   types_names = levels(types(object$configuration_list[[1]])),
  #                                   covariates_names = names(object$parameters$covariates),
  #                                   estimate_alpha = object$estimate_alpha,
  #                                   estimate_gamma = object$estimate_gamma)
  # y$se_numerical <- format_coefficient_vector(coefficient_vector = se_numerical,
  #                                             number_types = nlevels(types(object$configuration_list[[1]])),
  #                                             types_names = levels(types(object$configuration_list[[1]])),
  #                                             covariates_names = names(object$parameters$covariates),
  #                                             estimate_alpha = object$estimate_alpha,
  #                                             estimate_gamma = object$estimate_gamma)
  # y$lo <- format_coefficient_vector(coefficient_vector = lo,
  #                                   number_types = nlevels(types(object$configuration_list[[1]])),
  #                                   types_names = levels(types(object$configuration_list[[1]])),
  #                                   covariates_names = names(object$parameters$covariates),
  #                                   estimate_alpha = object$estimate_alpha,
  #                                   estimate_gamma = object$estimate_gamma)
  # y$hi <- format_coefficient_vector(coefficient_vector = hi,
  #                                   number_types = nlevels(types(object$configuration_list[[1]])),
  #                                   types_names = levels(types(object$configuration_list[[1]])),
  #                                   covariates_names = names(object$parameters$covariates),
  #                                   estimate_alpha = object$estimate_alpha,
  #                                   estimate_gamma = object$estimate_gamma)
  # y$lo_numerical <- format_coefficient_vector(coefficient_vector = lo_numerical,
  #                                             number_types = nlevels(types(object$configuration_list[[1]])),
  #                                             types_names = levels(types(object$configuration_list[[1]])),
  #                                             covariates_names = names(object$parameters$covariates),
  #                                             estimate_alpha = object$estimate_alpha,
  #                                             estimate_gamma = object$estimate_gamma)
  # y$hi_numerical <- format_coefficient_vector(coefficient_vector = hi_numerical,
  #                                             number_types = nlevels(types(object$configuration_list[[1]])),
  #                                             types_names = levels(types(object$configuration_list[[1]])),
  #                                             covariates_names = names(object$parameters$covariates),
  #                                             estimate_alpha = object$estimate_alpha,
  #                                             estimate_gamma = object$estimate_gamma)
  # y
}

#' Print summary of a Fitted Gibbs point process.
#'
#' @param x A summary of a fitted model object.
#' @param ... Ignored.
#' @export
print.summary_gibbsm <- function(x, ...) {
  print(x$coefficients)
}
