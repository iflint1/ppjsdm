remove(list = ls())
library(ppjsdm)
library(spatstat)
set.seed(1)

window <- Rectangle_window(c(-1, 1), c(-1, 1))
nreplications <- 1
ntypes <- 2
beta0 <- rep(2, ntypes)
estimate_gamma <- TRUE
alpha <- cbind(c(-0.3, 0.1), c(0.1, -0.2))

covariates <- list(temperature = function(x, y) x, elevation = function(x, y) y)
beta <- cbind(rep(2, ntypes), rep(1, ntypes))
ndummy <- 500

model <- "square_bump"
medium_range_model <- "Geyer"
saturation <- 2

short_range <- matrix(0.1, ntypes, ntypes)
medium_range <- matrix(0.1, ntypes, ntypes)
if(estimate_gamma) {
  long_range <- matrix(0.15, ntypes, ntypes)
  is_in <- matrix(NA, nrow = nreplications, ncol = ntypes + ntypes * (ntypes + 1) + length(covariates) * ntypes)
  gamma <- cbind(c(-0.5, -0.2), c(-0.2, -0.3))
} else {
  long_range <- matrix(0.1, ntypes, ntypes)
  is_in <- matrix(NA, nrow = nreplications, ncol = ntypes + ntypes * (ntypes + 1) / 2 + length(covariates) * ntypes)
  gamma <- matrix(0., ntypes, ntypes)
}
estimates <- matrix(NA, nrow = nreplications, ncol = ncol(is_in))
true <- c(beta0, alpha[1, 1], alpha[1, 2], alpha[2, 2])
if(estimate_gamma) {
  true <- c(true, gamma[1, 1], gamma[1, 2], gamma[2, 2])
}

if(length(covariates) > 0) {
  true <- c(true, beta[1, 1], beta[2, 1], beta[1, 2], beta[2, 2])
}

samples <- ppjsdm::rgibbs(window = window,
                          beta0 = beta0,
                          alpha = alpha,
                          gamma = gamma,
                          short_range = short_range,
                          medium_range = medium_range,
                          long_range = long_range,
                          nsim = nreplications,
                          model = model,
                          medium_range_model = medium_range_model,
                          drop = FALSE,
                          saturation = saturation,
                          covariates = covariates,
                          beta = beta)


for(i in seq_len(nreplications)) {
  fit <- gibbsm(samples[[i]],
                window = window,
                print = FALSE,
                short_range = short_range,
                medium_range = medium_range,
                long_range = long_range,
                use_glmnet = FALSE,
                ndummy = ndummy,
                model = model,
                medium_range_model = medium_range_model,
                covariates = covariates,
                saturation = saturation)

  s <- summary(fit)$coefficients
  estimate <- s$coefficients

  lower <- s$CI95_lo
  upper <- s$CI95_hi
  is_in[i, ] <- true >= lower & true <= upper
  estimates[i, ] <- estimate
}

print(colMeans(is_in, na.rm = TRUE))
print(colMeans(estimates, na.rm = TRUE))
