#' Calculate the matrix G2 used to compute the Variance-Covariance Matrix for a Fitted Gibbs point process.
#' The function can be supplied multiple fits, in which case the aggregate estimator is computed.
#' This is an internal function used when you want to compute the vcov matrix in multiple steps (useful for large datasets).
#'
#' @param ... A sequence of fit objects.
#' @param list List of fits.
#' @param nthreads (optional) number of threads to use.
#' @param debug Display debug information?
#' @export
compute_G2 <- function(..., list, nthreads, debug = FALSE) {
  if(missing(list)) {
    fits <- base::list(...)
  } else {
    fits <- list
  }

  theta <- setNames(sapply(seq_len(length(fits[[1]]$coefficients_vector)), function(i) mean(sapply(fits, function(fit) fit$coefficients_vector[i]), na.rm = TRUE)), nm = names(fits[[1]]$coefficients_vector))

  if(missing(nthreads)) {
    nthreads <- NULL
  }
  G2 <- lapply(fits, function(fit) {
    if(is.null(nthreads)) {
      nt <- fit$nthreads
    } else {
      nt <- nthreads
    }

    tt <- tryCatch({
      regressors <- as.matrix(fit$data_list$regressors)
    }, error = function(e) e, warning = function(w) w)

    # If we don't have enough RAM to call Matrix::as.matrix.Matrix, use custom function instead.
    if(is(tt, "error")) {
      # Ref: https://programmerah.com/the-sparse-matrix-of-r-language-is-too-large-to-be-used-as-matrix-8856/
      as_matrix <- function(mat) {
        tmp <- matrix(data = 0L, nrow = mat@Dim[1], ncol = mat@Dim[2])

        row_pos <- mat@i + 1
        col_pos <- findInterval(seq(mat@x) - 1, mat@p[-1]) + 1
        val <- mat@x

        for(i in seq_along(val)) {
          tmp[row_pos[i], col_pos[i]] <- val[i]
        }

        row.names(tmp) <- mat@Dimnames[[1]]
        colnames(tmp) <- mat@Dimnames[[2]]
        tmp
      }
      regressors <- as_matrix(fit$data_list$regressors)
    }

    if(debug) {
      cat(paste0("Starting computation of G2.\n"))
      tm <- Sys.time()
    }
    G2 <- compute_G2_cpp(configuration = fit$configuration_list[[1]],
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
                         regressors = regressors,
                         type = fit$data_list$type,
                         estimate_alpha = fit$estimate_alpha,
                         estimate_gamma = fit$estimate_gamma,
                         nthreads = nt,
                         dummy_distribution = fit$dummy_distribution,
                         mark_range = fit$mark_range)
    if(debug) {
      cat("End of computation. ")
      print(Sys.time() - tm)
    }
    G2
  })

  G2 <- Reduce("+", G2) / length(G2) / length(G2)

  rownames(G2) <- names(theta)
  colnames(G2) <- names(theta)

  G2
}
