remove(list = ls())
library(ppjsdm)
library(spatstat)

window <- Rectangle_window(c(0, 1), c(0, 1))
nreplications <- 1000
beta0 <- 4
alpha <- -0.3
ndummy <- 500
estimate_alpha <- TRUE
model <- "square_bump"
saturation <- 2
beta0_CI <- matrix(NA, nrow = nreplications, ncol = 2)
beta0_estimate <- vector(mode = "numeric", length = nreplications)
if(estimate_alpha) {
  is_in <- matrix(NA, nrow = nreplications, ncol = 2)
  alpha_CI <- matrix(NA, nrow = nreplications, ncol = 2)
  alpha_estimate <- vector(mode = "numeric", length = nreplications)
  short_range <- 0.1
} else {
  is_in <- matrix(NA, nrow = nreplications, ncol = 1)
  short_range <- 0.
}

# samples <- ppjsdm::rppp(window = window,
#                         lambda = exp(beta0),
#                         nsim = nreplications,
#                         drop = FALSE)

samples <- ppjsdm::rgibbs(window = window,
                          beta0 = beta0,
                          alpha = alpha,
                          short_range = short_range,
                          nsim = nreplications,
                          model = model,
                          drop = FALSE,
                          saturation = saturation)

# samples <- spatstat::rStrauss(beta = exp(beta0),
#                               gamma = exp(2 * alpha),
#                               R = short_range, W = as.owin(window),
#                               nsim = nreplications,
#                               drop = FALSE)
samples <- lapply(samples, function(s) ppjsdm::Configuration(x = s$x, y = s$y))


for(i in seq_len(nreplications)) {
  fit <- gibbsm(samples[[i]],
                window = window,
                print = FALSE,
                short_range = short_range,
                use_glmnet = FALSE,
                ndummy = ndummy,
                model = model,
                saturation = saturation)

  vc <- vcov(fit)

  se <- sqrt(diag(vc))
  estimate <- fit$coefficients

  if(estimate_alpha) {
    lower <- c(estimate$beta0 - 1.96 * se[1], estimate$alpha - 1.96 * se[2])
    upper <- c(estimate$beta0 + 1.96 * se[1], estimate$alpha + 1.96 * se[2])
    beta0_CI[i, ] <- c(lower[1], upper[1])
    alpha_CI[i, ] <- c(lower[2], upper[2])
    beta0_estimate[i] <- estimate$beta0
    alpha_estimate[i] <- estimate$alpha
    is_in[i, ] <- c(beta0, alpha) >= lower & c(beta0, alpha) <= upper
  } else {
    lower <- estimate$beta0 - 1.96 * se[1]
    upper <- estimate$beta0 + 1.96 * se[1]
    beta0_CI[i, ] <- c(lower, upper)
    beta0_estimate[i] <- estimate$beta0
    is_in[i, ] <- beta0 >= lower & beta0 <= upper
  }
}

print(colMeans(is_in, na.rm = TRUE))
if(estimate_alpha) {
  print(mean(alpha_estimate))
}
print(mean(beta0_estimate))
