#' Get parameters from a fit
#'
#' @param fit Fit.
#' @export
get_parameters_from_fit <- function(fit) {
  coefs <- fit$coefficients
  log_lambda_indices <- grep("log_lambda", rownames(coefs), value = TRUE)
  nspecies <- length(log_lambda_indices)

  lambda <- exp(coefs[log_lambda_indices, 1])

  alpha <- matrix(NA, nspecies, nspecies)
  for(i in seq_len(nspecies)) {
    for(j in seq_len(nspecies)) {
      alpha[i, j] <- coefs[paste0("alpha_", min(i, j), "_", max(i, j)), 1]
    }
  }

  gamma <- matrix(NA, nspecies, nspecies)
  for(i in seq_len(nspecies)) {
    for(j in seq_len(nspecies)) {
      gamma[i, j] <- coefs[paste0("gamma_", min(i, j), "_", max(i, j)), 1]
    }
  }

  list(lambda = lambda, alpha = alpha, gamma = gamma)
}
