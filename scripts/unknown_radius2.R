remove(list = ls())
library(ppjsdm)
library(spatstat)
set.seed(1)

window <- Rectangle_window(c(0, 1), c(0, 1))
nreplications <- 100
ntypes <- 2
beta0 <- c(3, 2.5)
steps <- 100000
estimate_gamma <- FALSE
alpha <- cbind(c(0.4, -0.6), c(-0.6, 0.4))

covariates <- list(temperature = function(x, y) x)
beta <- cbind(c(1.5, 2))
ndummy <- 1000

model <- "Geyer"
medium_range_model <- "square_exponential"
saturation <- 2

short_range <- cbind(c(0.04, 0.05), c(0.05, 0.06))
medium_range <- cbind(c(0.07, 0.02), c(0.02, 0.06))
if(estimate_gamma) {
  long_range <- cbind(c(0.17, 0.06), c(0.06, 0.2))
  is_in <- matrix(NA, nrow = nreplications, ncol = ntypes + ntypes * (ntypes + 1) + length(covariates) * ntypes)
  gamma <- cbind(c(-0.6, -0.8), c(-0.8, -0.4))
} else {
  long_range <- medium_range
  is_in <- matrix(NA, nrow = nreplications, ncol = ntypes + ntypes * (ntypes + 1) / 2 + length(covariates) * ntypes)
  gamma <- matrix(0., ntypes, ntypes)
}
radius_range <- c(0, 0.1)
estimates <- matrix(NA, nrow = nreplications, ncol = ncol(is_in))
true <- c(beta0, alpha[1, 1], alpha[1, 2], alpha[2, 2])
if(estimate_gamma) {
  true <- c(true, gamma[1, 1], gamma[1, 2], gamma[2, 2])
}

if(length(covariates) > 0) {
  true <- c(true, beta[1, 1], beta[2, 1])
}

short_range_estimates <- vector(mode = "list", length = nreplications)

tm <- Sys.time()

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

Sys.time() - tm

tm <- Sys.time()

for(i in seq_len(nreplications)) {
  fit <- gibbsm(samples[[i]],
                window = window,
                print = FALSE,
                short_range = radius_range,
                medium_range = c(0, 0),
                long_range = c(0, 0),
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
  short_range_estimates[[i]] <- fit$coefficients$short_range
}

Sys.time() - tm

coverage_probabilities <- colMeans(is_in, na.rm = TRUE)
mean_estimates <- colMeans(estimates, na.rm = TRUE)
median_estimates <- sapply(seq_len(ncol(estimates)), function(n) median(estimates[, n]))
mode_estimates <- sapply(seq_len(ncol(estimates)), function(n) rethinking::chainmode(estimates[, n]))
short_range_estimate <- Reduce("+", short_range_estimates) / length(short_range_estimates)

cum_short_range_estimates <- lapply(seq_len(nreplications), function(n) Reduce("+", short_range_estimates[1:n]) / n)
error <- sapply(cum_short_range_estimates, function(e) norm(e - short_range, type = "O"))
x <- 1:nreplications
lo <- loess(error ~ x, span = 0.2)
plot(x, error)
lines(predict(lo), col = 'red', lwd = 2)

average_error <- sapply(cum_short_range_estimates, function(e) mean(e - short_range, type = "O"))
x <- 1:nreplications
lo <- loess(average_error ~ x, span = 0.2)
plot(x, average_error)
lines(predict(lo), col = 'red', lwd = 2)

print(short_range_estimate)

output_string <- ""
output_string <- paste0(output_string, "$\\beta_{1,0}$ & ",
                        sprintf("%.2e", beta0[1]), " & ",
                        sprintf("%.2e", mean_estimates[1]), " & ",
                        sprintf("%.2f", coverage_probabilities[1]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{2,0}$ & ",
                        sprintf("%.2e", beta0[2]), " & ",
                        sprintf("%.2e", mean_estimates[2]), " & ",
                        sprintf("%.2f", coverage_probabilities[2]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\alpha_{1,1}$ & ",
                        sprintf("%.2e", alpha[1, 1]), " & ",
                        sprintf("%.2e", mean_estimates[3]), " & ",
                        sprintf("%.2f", coverage_probabilities[3]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\alpha_{1,2}$ & ",
                        sprintf("%.2e", alpha[1, 2]), " & ",
                        sprintf("%.2e", mean_estimates[4]), " & ",
                        sprintf("%.2f", coverage_probabilities[4]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\alpha_{2,2}$ & ",
                        sprintf("%.2e", alpha[2, 2]), " & ",
                        sprintf("%.2e", mean_estimates[5]), " & ",
                        sprintf("%.2f", coverage_probabilities[5]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{1,1}$ & ",
                        sprintf("%.2e", beta[1, 1]), " & ",
                        sprintf("%.2e", mean_estimates[6]), " & ",
                        sprintf("%.2f", coverage_probabilities[6]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{2,1}$ & ",
                        sprintf("%.2e", beta[2, 1]), " & ",
                        sprintf("%.2e", mean_estimates[7]), " & ",
                        sprintf("%.2f", coverage_probabilities[7]), " \\\\", sep = "\n")



output_string <- gsub("e\\+00", "", output_string)
cat(output_string)
