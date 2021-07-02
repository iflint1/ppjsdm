#' Calculate Variance-Covariance Matrix for a Fitted Gibbs point process.
#'
#' @param object A fitted model object.
#' @param npoints Target number of points in the restricted window that the vcov matrix is computed on. Computation is slower for larger values, but the vcov matrix is then better approximated.
#' @param multiple_windows Compute A2 and A3 on a lot of small windows and which are then averaged out, or only on a single restricted window?
#' @param ... Ignored.
#' @export
vcov.gibbsm <- function(object, npoints = 1000, multiple_windows = TRUE, ...) {
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
  vc$G1 + vc$G2
}
