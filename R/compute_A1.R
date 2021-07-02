#' Calculate the matrix A1 used to compute the Variance-Covariance Matrix for a Fitted Gibbs point process.
#' The function can be supplied multiple fits, in which case the aggregate estimator is computed.
#' This is an internal function used when you want to compute the vcov matrix in multiple steps (useful for large datasets).
#'
#' @param ... A sequence of fit objects.
#' @param list List of fits.
#' @param nthreads (optional) number of threads to use.
#' @param debug Display debug information?
#' @export
compute_A1 <- function(..., list, nthreads, debug = FALSE) {
  if(missing(list)) {
    fits <- base::list(...)
  } else {
    fits <- list
  }

  theta <- setNames(sapply(seq_len(length(fits[[1]]$coefficients_vector)), function(i) mean(sapply(fits, function(fit) fit$coefficients_vector[i]), na.rm = TRUE)), nm = names(fits[[1]]$coefficients_vector))

  if(missing(nthreads)) {
    nthreads <- NULL
  }
  A1 <- lapply(fits, function(fit) {
    if(is.null(nthreads)) {
      nt <- fit$nthreads
    } else {
      nt <- nthreads
    }
    if(debug) {
      cat(paste0("Starting computation of A1.\n"))
      tm <- Sys.time()
    }
    A1 <- compute_A1_cpp(rho = exp(-fit$data_list$shift),
                         theta = theta,
                         regressors = as.matrix(fit$data_list$regressors),
                         type = fit$data_list$type,
                         nthreads = nt)
    if(debug) {
      cat("End of computation. ")
      print(Sys.time() - tm)
    }
    A1
  })

  Reduce("+", A1) / length(A1)
}
