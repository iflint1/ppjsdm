#' Calculate the matrix G2 used to compute the Variance-Covariance Matrix for a Fitted Gibbs point process.
#' The function can be supplied multiple fits, in which case the aggregate estimator is computed.
#' This function is used to compute the vcov matrix in multiple steps (useful for large datasets).
#'
#' @param ... A sequence of fit objects.
#' @param list List of fits.
#' @param nthreads (optional) number of threads to use.
#' @param debug Display debug information?
#' @param time_limit Time limit measured in `unit` that can be spent running this function.
#' @param unit Unit used to measure the time limit (hours, mins, secs, etc).
#'
#' @export
compute_G2 <- function(..., list, nthreads = NULL, debug = FALSE, time_limit = Inf, unit = "hours") {
  # Allow for either sequence of fits or list of fits, convert both to list
  if(missing(list)) {
    fits <- base::list(...)
  } else {
    fits <- list
  }

  # Make sure thetas are compatible
  theta1 <- fits[[1]]$coefficients_vector
  for(fit in fits) {
    if(length(theta1) != length(fit$coefficients_vector)) {
      stop("Thetas of the supplied fits do not have the same length.")
    }
    if(!identical(names(theta1), names(fit$coefficients_vector))) {
      stop("Thetas of the supplied fits do not have the same names.")
    }
  }

  # Compute the regression coefficient, averaged out over the fits
  average_theta <- sapply(names(theta1), function(nm) {
    mean(sapply(fits, function(fit) fit$coefficients_vector[nm]), na.rm = TRUE)
  })

  compute_G2_on_fit <- function(fit) {
    # Use either the fit-specific nthreads, or the user-supplied value
    if(is.null(nthreads)) {
      if(is.null(fit$nthreads)) {
        nt <- 1
      } else {
        nt <- fit$nthreads
      }
    } else {
      nt <- nthreads
    }

    # Try to convert the regression matrix to base::matrix, it might crash if insufficient RAM
    tt <- tryCatch({
      regressors <- as.matrix(fit$data_list$regressors)
    }, error = function(e) e, warning = function(w) w)

    # If we don't have enough RAM to call Matrix::as.matrix.Matrix, use custom function instead.
    if(is(tt, "error") | is(tt, "warning")) {
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
      if(!inherits(fit$data_list$regressors, "Matrix")) {
        stop("Error while converting regression matrix to base::matrix format.")
      } else {
        regressors <- as_matrix(fit$data_list$regressors)
      }
    }

    if(debug) {
      cat(paste0("Starting computation of G2.\n"))
      current_time <- Sys.time()
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
                         theta = average_theta,
                         regressors = regressors,
                         type = fit$data_list$type,
                         estimate_alpha = fit$estimate_alpha,
                         estimate_gamma = fit$estimate_gamma,
                         nthreads = nt,
                         dummy_distribution = fit$dummy_distribution,
                         mark_range = fit$mark_range)
    if(debug) {
      cat(paste0("End of computation. Elapsed time: ", base::format(Sys.time() - current_time), ".\n"))
    }
    G2
  }

  # If no time limit supplied, compute G2 for each of the fits
  # and if not, try to guesstimate how long next fit will take
  # and execute as many as possible
  G2 <- execute_until_time_limit(objects = fits, func = compute_G2_on_fit, time_limit = time_limit, unit = unit)

  # Average out the matrices G2 over the fits and set names
  G2 <- Reduce("+", G2) / length(G2) / length(fits)
  rownames(G2) <- colnames(G2) <- names(average_theta)

  G2
}
