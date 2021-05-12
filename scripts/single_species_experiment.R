  remove(list = ls())
library(ppjsdm)
library(spatstat)
set.seed(1)

window <- Rectangle_window(c(0, 1), c(0, 1))
nreplications <- 200
ntypes <- 6
beta0 <- rep(3, ntypes)
steps <- 10000
estimate_gamma <- TRUE
alpha <- matrix(runif(1:(ntypes * ntypes), -2, 2), ntypes, ntypes)
alpha[lower.tri(alpha)] = t(alpha)[lower.tri(alpha)]
print(alpha)
diag_alpha <- diag(alpha)
covariates <- list(temperature = function(x, y) x, elevation = function(x, y) y)
beta <- cbind(rep(2, ntypes), rep(1, ntypes))
ndummy <- 1000

model <- "square_bump"
medium_range_model <- "square_exponential"
saturation <- 2

short_range <- matrix(0.01, ntypes, ntypes)
medium_range <- matrix(0.05, ntypes, ntypes)
long_range <- matrix(0.05, ntypes, ntypes)
is_in <- lapply(1:ntypes, function(n) matrix(NA, nrow = nreplications, ncol = 2 + length(covariates)))
gamma <- matrix(0., ntypes, ntypes)

estimates <- lapply(1:ntypes, function(n) matrix(NA, nrow = nreplications, ncol = ncol(is_in[[1]])))
true <- c(beta0, alpha[1, 1], alpha[1, 2], alpha[2, 2])

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
                          beta = beta,
                          steps = steps)


for(i in seq_len(nreplications)) {
  for(j in seq_len(ntypes)) {
    fit <- gibbsm(samples[[i]][j],
                  window = window,
                  short_range = short_range,
                  medium_range = medium_range,
                  long_range = long_range,
                  fitting_package = 'glm',
                  ndummy = ndummy,
                  model = model,
                  medium_range_model = medium_range_model,
                  covariates = covariates,
                  saturation = saturation)

    s <- summary(fit)$coefficients
    estimate <- s$coefficients

    lower <- s$CI95_lo
    upper <- s$CI95_hi
    val <- c(beta0[j], diag_alpha[j], beta[j, 1], beta[j, 2])
    is_in[[j]][i, ] <- val >= lower & val <= upper
    estimates[[j]][i, ] <- estimate
  }
}

print(sapply(is_in, function(n) colMeans(n, na.rm = TRUE)))
print(sapply(estimates, function(n) colMeans(n, na.rm = TRUE)))
