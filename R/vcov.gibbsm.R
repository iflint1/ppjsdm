#' Calculate Variance-Covariance Matrix for a Fitted Gibbs point process.
#'
#' @param object A fitted model object.
#' @param ... Ignored.
#' @export
vcov.gibbsm <- function(object, ...) {
  if(object$used_regularization) {
    stop("Can only compute VCOV matrix when no regularization was used.")
  }
  if(length(object$configuration_list) != 1) {
    stop("Cannot compute VCOV matrix for a fit obtained on a list of configurations.")
  }

  vc <- compute_vcov(configuration = object$configuration_list[[1]],
                     covariates = object$parameters$covariates,
                     model = object$parameters$model,
                     medium_range_model = object$parameters$medium_range_model,
                     short_range = object$coefficients$short_range,
                     medium_range = object$coefficients$medium_range,
                     long_range = object$coefficients$long_range,
                     saturation = object$parameters$saturation,
                     rho = exp(-object$data_list$shift),
                     window_volume = area(object$window),
                     theta = object$coefficients_vector,
                     regressors = object$data_list$regressors,
                     data_list = object$data_list,
                     estimate_alpha = object$estimate_alpha,
                     estimate_gamma = object$estimate_gamma,
                     nthreads = object$nthreads)
  vc$G1 + vc$G2
}
