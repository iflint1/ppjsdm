#' Calculate the matrix S used to compute the Variance-Covariance Matrix for a Fitted Gibbs point process.
#' The function can be supplied multiple fits, in which case the aggregate estimator is computed.
#' This is an internal function used when you want to compute the vcov matrix in multiple steps (useful for large datasets).
#'
#' @param ... A sequence of fit objects.
#' @param list List of fits.
#' @param nthreads (optional) number of threads to use.
#' @param debug Display debug information?
#' @param time_limit Time limit in hours that can be spent running this function.
#' @export
compute_S <- function(..., list, nthreads, debug = FALSE, time_limit) {
  if(missing(list)) {
    fits <- base::list(...)
  } else {
    fits <- list
  }

  theta <- setNames(sapply(seq_len(length(fits[[1]]$coefficients_vector)), function(i) mean(sapply(fits, function(fit) fit$coefficients_vector[i]), na.rm = TRUE)), nm = names(fits[[1]]$coefficients_vector))

  if(missing(nthreads)) {
    nthreads <- NULL
  }
  compute_S_on_fit <- function(fit) {
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
      cat(paste0("Starting computation of S.\n"))
      tm <- Sys.time()
    }
    S <- compute_S_cpp(rho = exp(-fit$data_list$shift),
                       theta = theta,
                       regressors = regressors,
                       type = fit$data_list$type,
                       nthreads = nt)
    if(debug) {
      cat("End of computation. ")
      print(Sys.time() - tm)
    }
    S
  }

  if(missing(time_limit)) {
    S <- lapply(fits, compute_S_on_fit)
  } else {
    S <- vector(mode = 'list', length = length(fits))
    execution_times <- c()
    time_left <- time_limit
    for(i in seq_len(length(fits))) {
      if(length(execution_times) > 0) {
        if(length(execution_times) > 1) {
          max_time <- mean(execution_times) + 4 * sd(execution_times)
        } else {
          max_time <- 2 * mean(execution_times)
        }
      } else {
        max_time <- 0
      }
      if(max_time < time_left) {
        start <- Sys.time()
        S[[i]] <- compute_S_on_fit(fits[[i]])
        execution_time <- as.numeric(Sys.time() - start, unit = "hours")
        execution_times <- c(execution_times, execution_time)
        time_left <- time_left - execution_time
      }
    }
    S <- S[sapply(S, function(s) !is.null(s))]
  }

  S <- Reduce("+", S) / length(S)

  rownames(S) <- names(theta)
  colnames(S) <- names(theta)

  S
}
