#' Summary of a Fitted Gibbs point process.
#'
#' @param object A fitted model object.
#' @param npoints Target number of points in the restricted window that the vcov matrix is computed on. Computation is slower for larger values, but the vcov matrix is then better approximated.
#' @param multiple_windows Compute A2 and A3 on a lot of small windows and which are then averaged out, or only on a single restricted window?
#' @param ... Ignored.
#' @importFrom stats pnorm qnorm
#' @export
summary.gibbsm <- function(object, npoints = 1000, multiple_windows = TRUE, ...) {
  y <- list()
  class(y) <- "summary_gibbsm"

  if(object$used_regularization) {
    warning("Computing the Variance-Covariance matrix of a regularised fit.")
  }
  if(length(object$configuration_list) != 1) {
    stop("Cannot compute VCOV matrix for a fit obtained on a list of configurations.")
  }

  vc <- compute_vcov(configuration = object$configuration_list[[1]],
                     dummy = object$data_list$dummy,
                     window = object$window,
                     covariates = object$parameters$covariates,
                     model = object$parameters$model,
                     medium_range_model = object$parameters$medium_range_model,
                     short_range = object$coefficients$short_range,
                     medium_range = object$coefficients$medium_range,
                     long_range = object$coefficients$long_range,
                     saturation = object$parameters$saturation,
                     rho = exp(-object$data_list$shift),
                     theta = object$coefficients_vector,
                     regressors = as.matrix(object$data_list$regressors),
                     data_list = object$data_list,
                     estimate_alpha = object$estimate_alpha,
                     estimate_gamma = object$estimate_gamma,
                     debug = object$debug,
                     nthreads = object$nthreads,
                     npoints = npoints,
                     multiple_windows = multiple_windows,
                     dummy_distribution = object$dummy_distribution,
                     mark_range = object$mark_range)

  se_numerical <- sqrt(diag(vc$G2))
  se <- sqrt(diag(vc$G1 + vc$G2))

  se_numerical_proportion <- se_numerical / se
  coefficients <- object$coefficients_vector
  one_ninetysix <- qnorm(0.975)
  lo <- coefficients - one_ninetysix * se
  hi <- coefficients + one_ninetysix * se
  lo_numerical <- coefficients - one_ninetysix * se_numerical
  hi_numerical <- coefficients + one_ninetysix * se_numerical
  zval <- coefficients / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  psig <- cut(pval,
              breaks = c(0,0.001, 0.01, 0.05, 1),
              labels = c("***", "**", "*", "  "),
              include.lowest = TRUE)
  y$G1 <- vc$G1
  y$G2 <- vc$G2
  y$coefficients <- data.frame(coefficients = coefficients,
                               se = se,
                               CI95_lo = lo,
                               CI95_hi = hi,
                               Ztest = psig,
                               Pval = pval,
                               Zval = zval,
                               se_numerical_proportion = se_numerical_proportion)

  y$se <- format_coefficient_vector(coefficient_vector = se,
                                    number_types = nlevels(types(object$configuration_list[[1]])),
                                    types_names = levels(types(object$configuration_list[[1]])),
                                    covariates_names = names(object$parameters$covariates),
                                    estimate_alpha = object$estimate_alpha,
                                    estimate_gamma = object$estimate_gamma)
  y$se_numerical <- format_coefficient_vector(coefficient_vector = se_numerical,
                                              number_types = nlevels(types(object$configuration_list[[1]])),
                                              types_names = levels(types(object$configuration_list[[1]])),
                                              covariates_names = names(object$parameters$covariates),
                                              estimate_alpha = object$estimate_alpha,
                                              estimate_gamma = object$estimate_gamma)
  y$lo <- format_coefficient_vector(coefficient_vector = lo,
                                    number_types = nlevels(types(object$configuration_list[[1]])),
                                    types_names = levels(types(object$configuration_list[[1]])),
                                    covariates_names = names(object$parameters$covariates),
                                    estimate_alpha = object$estimate_alpha,
                                    estimate_gamma = object$estimate_gamma)
  y$hi <- format_coefficient_vector(coefficient_vector = hi,
                                    number_types = nlevels(types(object$configuration_list[[1]])),
                                    types_names = levels(types(object$configuration_list[[1]])),
                                    covariates_names = names(object$parameters$covariates),
                                    estimate_alpha = object$estimate_alpha,
                                    estimate_gamma = object$estimate_gamma)
  y$lo_numerical <- format_coefficient_vector(coefficient_vector = lo_numerical,
                                              number_types = nlevels(types(object$configuration_list[[1]])),
                                              types_names = levels(types(object$configuration_list[[1]])),
                                              covariates_names = names(object$parameters$covariates),
                                              estimate_alpha = object$estimate_alpha,
                                              estimate_gamma = object$estimate_gamma)
  y$hi_numerical <- format_coefficient_vector(coefficient_vector = hi_numerical,
                                              number_types = nlevels(types(object$configuration_list[[1]])),
                                              types_names = levels(types(object$configuration_list[[1]])),
                                              covariates_names = names(object$parameters$covariates),
                                              estimate_alpha = object$estimate_alpha,
                                              estimate_gamma = object$estimate_gamma)
  y
}

#' Print summary of a Fitted Gibbs point process.
#'
#' @param x A summary of a fitted model object.
#' @param ... Ignored.
#' @export
print.summary_gibbsm <- function(x, ...) {
  print(x$coefficients)
}
