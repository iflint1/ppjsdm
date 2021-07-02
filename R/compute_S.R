#' Calculate the matrix S used to compute the Variance-Covariance Matrix for a Fitted Gibbs point process.
#' The function can be supplied multiple fits, in which case the aggregate estimator is computed.
#' This is an internal function used when you want to compute the vcov matrix in multiple steps (useful for large datasets).
#'
#' @param ... A sequence of fit objects.
#' @param list List of fits.
#' @export
compute_S <- function(..., list) {
  if(missing(list)) {
    fits <- list(...)
  } else {
    fits <- list
  }

  theta <- setNames(sapply(seq_len(length(fits[[1]]$coefficients_vector)), function(i) mean(sapply(fits, function(fit) fit$coefficients_vector[i]), na.rm = TRUE)), nm = names(fits[[1]]$coefficients_vector))

  S <- lapply(fits, function(fit) {
    compute_S_cpp(rho = exp(-fit$data_list$shift),
                  theta = theta,
                  regressors = as.matrix(fit$data_list$regressors),
                  type = fit$data_list$type,
                  nthreads = fit$nthreads)
  })

  Reduce("+", S) / length(S)
}
