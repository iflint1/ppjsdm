#' Calculate the matrix G2 used to compute the Variance-Covariance Matrix for a Fitted Gibbs point process.
#' The function can be supplied multiple fits, in which case the aggregate estimator is computed.
#' This is an internal function used when you want to compute the vcov matrix in multiple steps (useful for large datasets).
#'
#' @param ... A list of fit objects.
#' @export
compute_G2 <- function(..., npoints = 1000, multiple_windows = TRUE) {
  fits <- list(...)

  theta <- setNames(sapply(seq_len(length(fits[[1]]$coefficients_vector)), function(i) mean(sapply(fits, function(fit) fit$coefficients_vector[i]), na.rm = TRUE)), nm = names(fits[[1]]$coefficients_vector))

  G2 <- lapply(fits, function(fit) {
    compute_G2_cpp(configuration = fit$configuration_list[[1]],
                   dummy = fit$data_list$dummy,
                   window = fit$window,
                   covariates = fit$parameters$covariates,
                   model = fit$parameters$model,
                   medium_range_model = fit$parameters$medium_range_model,
                   short_range = fit$coefficients$short_range,
                   medium_range = fit$coefficients$medium_range,
                   long_range = fit$coefficients$long_range,
                   saturation = fit$parameters$saturation,
                   rho = exp(-fit$data_list$shift),
                   theta = theta,
                   regressors = as.matrix(fit$data_list$regressors),
                   type = fit$data_list$type,
                   estimate_alpha = fit$estimate_alpha,
                   estimate_gamma = fit$estimate_gamma,
                   nthreads = fit$nthreads,
                   dummy_distribution = fit$dummy_distribution,
                   mark_range = fit$mark_range)
  })

  Reduce("+", G2) / length(G2) / length(G2)
}
