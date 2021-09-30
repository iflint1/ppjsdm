#' Calculate the matrix A_2 + A3 used to compute the Variance-Covariance Matrix for a Fitted Gibbs point process.
#' The function can be supplied multiple fits, in which case the aggregate estimator is computed.
#' This is an internal function used when you want to compute the vcov matrix in multiple steps (useful for large datasets).
#'
#' @param ... A sequence of fit objects.
#' @param list List of fits.
#' @param nthreads (optional) number of threads to use.
#' @param debug Display debug information?
#' @param npoints Target number of points in the restricted window that the vcov matrix is computed on. Computation is slower for larger values, but the vcov matrix is then better approximated.
#' @param multiple_windows Compute A2 and A3 on a lot of small windows and which are then averaged out, or only on a single restricted window?
#' @param time_limit Time limit in hours that can be spent running this function.
#' @importFrom stats sd
#' @export
compute_A2_plus_A3 <- function(..., list, nthreads, debug = FALSE, npoints = 1000, multiple_windows = TRUE, time_limit) {
  if(missing(list)) {
    fits <- base::list(...)
  } else {
    fits <- list
  }

  theta <- setNames(sapply(seq_len(length(fits[[1]]$coefficients_vector)), function(i) mean(sapply(fits, function(fit) fit$coefficients_vector[i]), na.rm = TRUE)), nm = names(fits[[1]]$coefficients_vector))

  if(missing(nthreads)) {
    nthreads <- fits[[1]]$nthreads
  }

  tt <- tryCatch({
    regressors <- as.matrix(fits[[1]]$data_list$regressors)
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
    regressors <- as_matrix(fits[[1]]$data_list$regressors)
  }

  if(debug) {
    cat(paste0("Starting computation of A2 + A3.\n"))
    tm <- Sys.time()
  }
  if(missing(time_limit)) {
    A2_plus_A3 <- compute_A2_plus_A3_cpp(configuration = fits[[1]]$configuration_list[[1]],
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
                                         regressors = regressors,
                                         data_list = fits[[1]]$data_list,
                                         estimate_alpha = fits[[1]]$estimate_alpha,
                                         estimate_gamma = fits[[1]]$estimate_gamma,
                                         nthreads = nthreads,
                                         npoints = npoints,
                                         multiple_windows = multiple_windows,
                                         mark_range = fits[[1]]$mark_range,
                                         debug = debug,
                                         max_executions = 1e5)
  } else {
    start <- Sys.time()
    A2_plus_A3 <- compute_A2_plus_A3_cpp(configuration = fits[[1]]$configuration_list[[1]],
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
                                         regressors = regressors,
                                         data_list = fits[[1]]$data_list,
                                         estimate_alpha = fits[[1]]$estimate_alpha,
                                         estimate_gamma = fits[[1]]$estimate_gamma,
                                         nthreads = nthreads,
                                         npoints = npoints,
                                         multiple_windows = multiple_windows,
                                         mark_range = fits[[1]]$mark_range,
                                         debug = debug,
                                         max_executions = 1)
    single_execution_time <- as.numeric(Sys.time() - start, unit = "hours")
    max_executions <- max(0, floor(time_limit / single_execution_time - 1.))
    if(debug) {
      cat(paste0("Going to divide the window into a maximum of ", max_executions, " subregions to approximate A2 + A3.\n"))
    }
    if(multiple_windows & max_executions > 1) {
      A2_plus_A3 <- compute_A2_plus_A3_cpp(configuration = fits[[1]]$configuration_list[[1]],
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
                                           regressors = regressors,
                                           data_list = fits[[1]]$data_list,
                                           estimate_alpha = fits[[1]]$estimate_alpha,
                                           estimate_gamma = fits[[1]]$estimate_gamma,
                                           nthreads = nthreads,
                                           npoints = npoints,
                                           multiple_windows = multiple_windows,
                                           mark_range = fits[[1]]$mark_range,
                                           debug = debug,
                                           max_executions = max_executions)
    }
  }
  if(debug) {
    cat("End of computation. ")
    print(Sys.time() - tm)
  }

  rownames(A2_plus_A3) <- names(theta)
  colnames(A2_plus_A3) <- names(theta)
  A2_plus_A3
}
