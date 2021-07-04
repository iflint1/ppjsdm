#' Calculate the matrix S used to compute the Variance-Covariance Matrix for a Fitted Gibbs point process.
#' The function can be supplied multiple fits, in which case the aggregate estimator is computed.
#' This is an internal function used when you want to compute the vcov matrix in multiple steps (useful for large datasets).
#'
#' @param ... A sequence of fit objects.
#' @param list List of fits.
#' @param nthreads (optional) number of threads to use.
#' @param debug Display debug information?
#' @export
compute_S <- function(..., list, nthreads, debug = FALSE) {
  if(missing(list)) {
    fits <- base::list(...)
  } else {
    fits <- list
  }

  theta <- setNames(sapply(seq_len(length(fits[[1]]$coefficients_vector)), function(i) mean(sapply(fits, function(fit) fit$coefficients_vector[i]), na.rm = TRUE)), nm = names(fits[[1]]$coefficients_vector))

  # Get rid of everything superfluous to facilitate working with large matrices
  fits <- lapply(fits, function(fit) {
    base::list(rho = exp(-fit$data_list$shift),
               regressors = fit$data_list$regressors,
               type = fit$data_list$type,
               nthreads = fit$nthreads)
  })
  gc()

  if(missing(nthreads)) {
    nthreads <- NULL
  }

  S <- vector(mode = 'list', length = length(fits))
  for(i in seq_len(length(fits))) {
    gc()
    if(is.null(nthreads)) {
      nt <- fits[[i]]$nthreads
    } else {
      nt <- nthreads
    }
    if(debug) {
      cat(paste0("Starting computation of S.\n"))
      tm <- Sys.time()
    }
    S[[i]] <- compute_S_cpp(rho = fits[[i]]$rho,
                            theta = theta,
                            regressors = as.matrix(fits[[i]]$regressors),
                            type = fits[[i]]$type,
                            nthreads = nt)
    if(debug) {
      cat("End of computation. ")
      print(Sys.time() - tm)
    }
  }

  Reduce("+", S) / length(S)
}
