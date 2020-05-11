remove(list = ls())
library(spatstat)

window <- owin()
nreplications <- 10
log_beta <- 4
log_gamma <- -0.3
R = 0.1
log_beta_CI <- matrix(NA, nrow = nreplications, ncol = 2)
log_beta_estimate <- vector(mode = "numeric", length = nreplications)
is_in <- matrix(NA, nrow = nreplications, ncol = 2)
log_gamma_CI <- matrix(NA, nrow = nreplications, ncol = 2)
log_gamma_estimate <- vector(mode = "numeric", length = nreplications)

samples <- spatstat::rStrauss(beta = exp(log_beta), gamma = exp(log_gamma), R = R, W = window, nsim = nreplications, drop = FALSE)

for(i in seq_len(nreplications)) {
  Y <- samples[[i]]
  fit <- spatstat::ppm(Y ~ 1, method = "logi", interaction = Strauss(R), correction = "none")

  vc <- vcov(fit)
  se <- sqrt(diag(vc))
  estimate <- coefficients(fit)

  lower <- c(estimate[1] - 1.96 * se[1], estimate[2] - 1.96 * se[2])
  upper <- c(estimate[1] + 1.96 * se[1], estimate[2] + 1.96 * se[2])
  log_beta_CI[i, ] <- c(lower[1], upper[1])
  log_gamma_CI[i, ] <- c(lower[2], upper[2])
  log_beta_estimate[i] <- estimate[1]
  log_gamma_estimate[i] <- estimate[2]
  is_in[i, ] <- c(log_beta, log_gamma) >= lower & c(log_beta, log_gamma) <= upper
}

print(colMeans(is_in, na.rm = TRUE))
print(mean(log_gamma_estimate))
print(mean(log_beta_estimate))
