remove(list = ls())
library(ppjsdm)
library(spatstat)
set.seed(1)

window <- Rectangle_window(c(-1, 1), c(-1, 1))
nreplications <- 100
ntypes <- 2
beta0 <- c(3.5, 3)
steps <- 0
estimate_gamma <- TRUE
alpha <- cbind(c(-0.2, 0), c(0, -0.4))

covariates <- list(temperature = function(x, y) x)
beta <- cbind(c(1.5, 2))
ndummy <- 1000

model <- "square_bump"
medium_range_model <- "square_exponential"
saturation <- 2

short_range <- matrix(0.05, ntypes, ntypes)
medium_range <- matrix(0.05, ntypes, ntypes)
if(estimate_gamma) {
  long_range <- matrix(0.1, ntypes, ntypes)
  is_in <- matrix(NA, nrow = nreplications, ncol = ntypes + 2 * ntypes + length(covariates) * ntypes)
  gamma <- cbind(c(-0.6, 0), c(0, -0.1))
} else {
  long_range <- matrix(0.05, ntypes, ntypes)
  is_in <- matrix(NA, nrow = nreplications, ncol = ntypes + ntypes + length(covariates) * ntypes)
  gamma <- matrix(0., ntypes, ntypes)
}
estimates <- matrix(NA, nrow = nreplications, ncol = ncol(is_in))
true <- c(beta0, alpha[1, 1], alpha[2, 2])
if(estimate_gamma) {
  true <- c(true, gamma[1, 1], gamma[2, 2])
}

if(length(covariates) > 0) {
  true <- c(true, beta[1, 1], beta[2, 1])
}

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
  fit1 <- gibbsm(samples[[i]][1],
                 window = window,
                 short_range = short_range,
                 medium_range = medium_range,
                 long_range = long_range,
                 use_glmnet = FALSE,
                 ndummy = ndummy,
                 model = model,
                 medium_range_model = medium_range_model,
                 covariates = covariates,
                 saturation = saturation)
  fit2 <- gibbsm(samples[[i]][2],
                 window = window,
                 short_range = short_range,
                 medium_range = medium_range,
                 long_range = long_range,
                 use_glmnet = FALSE,
                 ndummy = ndummy,
                 model = model,
                 medium_range_model = medium_range_model,
                 covariates = covariates,
                 saturation = saturation)

  s1 <- summary(fit1)$coefficients
  estimate1 <- s1$coefficients

  s2 <- summary(fit2)$coefficients
  estimate2 <- s2$coefficients

  lower1 <- s1$CI95_lo
  upper1 <- s1$CI95_hi

  lower2 <- s2$CI95_lo
  upper2 <- s2$CI95_hi

  lower <- c(lower1[1], lower2[1], lower1[2], lower2[2], lower1[3], lower2[3], lower1[4], lower2[4])
  upper <- c(upper1[1], upper2[1], upper1[2], upper2[2], upper1[3], upper2[3], upper1[4], upper2[4])
  estimate <- c(estimate1[1], estimate2[1], estimate1[2], estimate2[2], estimate1[3], estimate2[3], estimate1[4], estimate2[4])

  is_in[i, ] <- true >= lower & true <= upper
  estimates[i, ] <- estimate
}

Sys.time() - tm

coverage_probabilities <- colMeans(is_in, na.rm = TRUE)
mean_estimates <- colMeans(estimates, na.rm = TRUE)
median_estimates <- sapply(seq_len(ncol(estimates)), function(n) median(estimates[, n]))
mode_estimates <- sapply(seq_len(ncol(estimates)), function(n) rethinking::chainmode(estimates[, n]))

output_string <- ""
output_string <- paste0(output_string, "$\\beta_{1,0}$ & ",
                        sprintf("%.2e", beta0[1]), " & ",
                        sprintf("%.2e", mean_estimates[1]), " & ",
                        sprintf("%.2e", median_estimates[1]), " & ",
                        sprintf("%.2e", mode_estimates[1]), " & ",
                        sprintf("%.2f", coverage_probabilities[1]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{2,0}$ & ",
                        sprintf("%.2e", beta0[2]), " & ",
                        sprintf("%.2e", mean_estimates[2]), " & ",
                        sprintf("%.2e", median_estimates[2]), " & ",
                        sprintf("%.2e", mode_estimates[2]), " & ",
                        sprintf("%.2f", coverage_probabilities[2]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\alpha_{1,1}$ & ",
                        sprintf("%.2e", alpha[1, 1]), " & ",
                        sprintf("%.2e", mean_estimates[3]), " & ",
                        sprintf("%.2e", median_estimates[3]), " & ",
                        sprintf("%.2e", mode_estimates[3]), " & ",
                        sprintf("%.2f", coverage_probabilities[3]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\alpha_{1,2}$ & ",
                        sprintf("%.2e", alpha[1, 2]), " & ",
                        sprintf("%.2e", mean_estimates[4]), " & ",
                        sprintf("%.2e", median_estimates[4]), " & ",
                        sprintf("%.2e", mode_estimates[4]), " & ",
                        sprintf("%.2f", coverage_probabilities[4]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\alpha_{2,2}$ & ",
                        sprintf("%.2e", alpha[2, 2]), " & ",
                        sprintf("%.2e", mean_estimates[5]), " & ",
                        sprintf("%.2e", median_estimates[5]), " & ",
                        sprintf("%.2e", mode_estimates[5]), " & ",
                        sprintf("%.2f", coverage_probabilities[5]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\gamma_{1,1}$ & ",
                        sprintf("%.2e", gamma[1, 1]), " & ",
                        sprintf("%.2e", mean_estimates[6]), " & ",
                        sprintf("%.2e", median_estimates[6]), " & ",
                        sprintf("%.2e", mode_estimates[6]), " & ",
                        sprintf("%.2f", coverage_probabilities[6]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\gamma_{1,2}$ & ",
                        sprintf("%.2e", gamma[1, 2]), " & ",
                        sprintf("%.2e", mean_estimates[7]), " & ",
                        sprintf("%.2e", median_estimates[7]), " & ",
                        sprintf("%.2e", mode_estimates[7]), " & ",
                        sprintf("%.2f", coverage_probabilities[7]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\gamma_{2,2}$ & ",
                        sprintf("%.2e", gamma[2, 2]), " & ",
                        sprintf("%.2e", mean_estimates[8]), " & ",
                        sprintf("%.2e", median_estimates[8]), " & ",
                        sprintf("%.2e", mode_estimates[8]), " & ",
                        sprintf("%.2f", coverage_probabilities[8]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{1,1}$ & ",
                        sprintf("%.2e", beta[1, 1]), " & ",
                        sprintf("%.2e", mean_estimates[9]), " & ",
                        sprintf("%.2e", median_estimates[9]), " & ",
                        sprintf("%.2e", mode_estimates[9]), " & ",
                        sprintf("%.2f", coverage_probabilities[9]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{2,1}$ & ",
                        sprintf("%.2e", beta[2, 1]), " & ",
                        sprintf("%.2e", mean_estimates[10]), " & ",
                        sprintf("%.2e", median_estimates[10]), " & ",
                        sprintf("%.2e", mode_estimates[10]), " & ",
                        sprintf("%.2f", coverage_probabilities[10]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{1,2}$ & ",
                        sprintf("%.2e", beta[1, 2]), " & ",
                        sprintf("%.2e", mean_estimates[11]), " & ",
                        sprintf("%.2e", median_estimates[11]), " & ",
                        sprintf("%.2e", mode_estimates[11]), " & ",
                        sprintf("%.2f", coverage_probabilities[11]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{2,2}$ & ",
                        sprintf("%.2e", beta[2, 2]), " & ",
                        sprintf("%.2e", mean_estimates[12]), " & ",
                        sprintf("%.2e", median_estimates[12]), " & ",
                        sprintf("%.2e", mode_estimates[12]), " & ",
                        sprintf("%.2f", coverage_probabilities[12]), " \\\\", sep = "\n")



output_string <- gsub("e\\+00", "", output_string)
cat(output_string)
