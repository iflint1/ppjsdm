remove(list = ls())
library(ppjsdm)

window <- Rectangle_window(c(0, 1), c(0, 1))
nreplications <- 1000
log_lambda <- 3
alpha <- 0
estimate_alpha <- FALSE
log_lambda_CI <- matrix(NA, nrow = nreplications, ncol = 2)
log_lambda_estimate <- vector(mode = "numeric", length = nreplications)
samples <- vector(mode = "list", length = nreplications)
if(estimate_alpha) {
  is_in <- matrix(NA, nrow = nreplications, ncol = 2)
  alpha_CI <- matrix(NA, nrow = nreplications, ncol = 2)
  alpha_estimate <- vector(mode = "numeric", length = nreplications)
  short_range <- 0.2
} else {
  is_in <- matrix(NA, nrow = nreplications, ncol = 1)
  short_range <- 0.
}

for(i in seq_len(nreplications)) {
  samples[[i]] <- ppjsdm::rgibbs(window = window, lambda = exp(log_lambda), alpha = alpha, short_range = short_range)
  fit <- gibbsm(samples[[i]], window = window, print = FALSE, short_range = short_range, use_glmnet = FALSE)

  vc <- vcov(fit)
  se <- sqrt(diag(vc))
  estimate <- fit$coefficients

  if(estimate_alpha) {
    lower <- c(log(estimate$lambda) - 1.96 * se[1], estimate$alpha - 1.96 * se[2])
    upper <- c(log(estimate$lambda) + 1.96 * se[1], estimate$alpha + 1.96 * se[2])
    log_lambda_CI[i, ] <- c(lower[1], upper[1])
    alpha_CI[i, ] <- c(lower[2], upper[2])
    log_lambda_estimate[i] <- log(estimate$lambda)
    alpha_estimate[i] <- estimate$alpha
    is_in[i, ] <- c(log_lambda, alpha) >= lower & c(log_lambda, alpha) <= upper
  } else {
    lower <- log(estimate$lambda) - 1.96 * se[1]
    upper <- log(estimate$lambda) + 1.96 * se[1]
    log_lambda_CI[i, ] <- c(lower, upper)
    log_lambda_estimate[i] <- log(estimate$lambda)
    is_in[i, ] <- log_lambda >= lower & log_lambda <= upper
  }
}

print(colMeans(is_in, na.rm = TRUE))
if(estimate_alpha) {
  print(mean(alpha_estimate))
}
print(mean(log_lambda_estimate))
