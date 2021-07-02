#' Calculate the matrix A_2 + A3 used to compute the Variance-Covariance Matrix for a Fitted Gibbs point process.
#' The function can be supplied multiple fits, in which case the aggregate estimator is computed.
#' This is an internal function used when you want to compute the vcov matrix in multiple steps (useful for large datasets).
#'
#' @param ... A sequence of fit objects.
#' @param list List of fits.
#' @param npoints Target number of points in the restricted window that the vcov matrix is computed on. Computation is slower for larger values, but the vcov matrix is then better approximated.
#' @param multiple_windows Compute A2 and A3 on a lot of small windows and which are then averaged out, or only on a single restricted window?
#' @export
compute_A2_plus_A3 <- function(..., npoints = 1000, multiple_windows = TRUE) {
  if(missing(list)) {
    fits <- list(...)
  } else {
    fits <- list
  }

  theta <- setNames(sapply(seq_len(length(fits[[1]]$coefficients_vector)), function(i) mean(sapply(fits, function(fit) fit$coefficients_vector[i]), na.rm = TRUE)), nm = names(fits[[1]]$coefficients_vector))

  compute_A2_plus_A3_cpp(configuration = fits[[1]]$configuration_list[[1]],
                         window = fits[[1]]$window,
                         covariates = fits[[1]]$parameters$covariates,
                         model = fits[[1]]$parameters$model,
                         medium_range_model = fits[[1]]$parameters$medium_range_model,
                         short_range = fits[[1]]$coefficients$short_range,
                         medium_range = fits[[1]]$coefficients$medium_range,
                         long_range = fits[[1]]$coefficients$long_range,
                         saturation = fits[[1]]$parameters$saturation,
                         rho = exp(-fits[[1]]$data_list$shift),
                         theta = theta,
                         regressors = as.matrix(fits[[1]]$data_list$regressors),
                         data_list = fits[[1]]$data_list,
                         estimate_alpha = fits[[1]]$estimate_alpha,
                         estimate_gamma = fits[[1]]$estimate_gamma,
                         nthreads = fits[[1]]$nthreads,
                         npoints = npoints,
                         multiple_windows = multiple_windows,
                         mark_range = fits[[1]]$mark_range)
}
