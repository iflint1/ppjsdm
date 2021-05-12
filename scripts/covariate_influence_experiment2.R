remove(list = ls())
library(ppjsdm)
library(spatstat)
set.seed(1)

window <- Rectangle_window(c(0, 1), c(0, 1))
nreplications <- 10
ntypes <- 2
steps <- 100000
beta0 <- c(3.8, 3.8)
alpha <- cbind(c(0, 1.5), c(1.5, 0))

covariates <- list(temperature = function(x, y) x)
beta <- matrix(c(2, 1), nrow = ntypes, ncol = 1)
ndummy <- 1000

model <- "square_bump"
medium_range_model <- "square_exponential"
saturation <- 2

short_range <- matrix(0.01, ntypes, ntypes)
medium_range <- matrix(0.05, ntypes, ntypes)
long_range <- matrix(0.05, ntypes, ntypes)
is_in <- matrix(NA, nrow = nreplications, ncol = ntypes + ntypes * (ntypes + 1) / 2 + length(covariates) * ntypes)
is_in_1 <- matrix(NA, nrow = nreplications, ncol = 2 + length(covariates))
is_in_2 <- matrix(NA, nrow = nreplications, ncol = ncol(is_in_1))
estimates <- matrix(NA, nrow = nreplications, ncol = ncol(is_in))
estimates_1 <- matrix(NA, nrow = nreplications, ncol = ncol(is_in_1))
estimates_2 <- matrix(NA, nrow = nreplications, ncol = ncol(is_in_2))
gamma <- matrix(0., ntypes, ntypes)


true <- c(beta0, alpha[1, 1], alpha[1, 2], alpha[2, 2])
true_1 <- c(beta0[1], alpha[1, 1])
true_2 <- c(beta0[2], alpha[2, 2])

if(length(covariates) > 0) {
  true <- c(true, beta[1, 1], beta[2, 1])
  true_1 <- c(true_1, beta[1, 1])
  true_2 <- c(true_2, beta[2, 1])
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
  print(i)
  fit <- gibbsm(samples[[i]],
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
  fit1 <- gibbsm(samples[[i]][1],
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
  fit2 <- gibbsm(samples[[i]][2],
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

  s1 <- summary(fit1)$coefficients
  estimate1 <- s1$coefficients

  s2 <- summary(fit2)$coefficients
  estimate2 <- s2$coefficients

  lower <- s$CI95_lo
  upper <- s$CI95_hi
  is_in[i, ] <- true >= lower & true <= upper
  estimates[i, ] <- estimate

  lower1 <- s1$CI95_lo
  upper1 <- s1$CI95_hi
  is_in_1[i, ] <- true_1 >= lower1 & true_1 <= upper1
  estimates_1[i, ] <- estimate1

  lower2 <- s2$CI95_lo
  upper2 <- s2$CI95_hi
  is_in_2[i, ] <- true_2 >= lower2 & true_2 <= upper2
  estimates_2[i, ] <- estimate2
}

Sys.time() - tm

RMSE <- function(estimate, true){
  sqrt(mean((estimate - true)^2))
}

coverage_probabilities <- colMeans(is_in, na.rm = TRUE)
mean_estimates <- colMeans(estimates, na.rm = TRUE)
min_estimates <- apply(estimates, 2, min)
max_estimates <- apply(estimates, 2, max)
rmse_estimates <- sapply(seq_len(ncol(estimates)), function(n) RMSE(estimates[, n], true[n]))

output_string <- ""
output_string <- paste0(output_string, "$\\beta_{1,0}$ & ",
                        sprintf("%.2e", beta0[1]), " & ",
                        sprintf("%.2e", mean_estimates[1]), " & ",
                        sprintf("%.2e", min_estimates[1]), " & ",
                        sprintf("%.2e", max_estimates[1]), " & ",
                        sprintf("%.2f", coverage_probabilities[1]), " & ",
                        sprintf("%.2e", rmse_estimates[1]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{2,0}$ & ",
                        sprintf("%.2e", beta0[2]), " & ",
                        sprintf("%.2e", mean_estimates[2]), " & ",
                        sprintf("%.2e", min_estimates[2]), " & ",
                        sprintf("%.2e", max_estimates[2]), " & ",
                        sprintf("%.2f", coverage_probabilities[2]), " & ",
                        sprintf("%.2e", rmse_estimates[2]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\alpha_{1,1}$ & ",
                        sprintf("%.2e", alpha[1, 1]), " & ",
                        sprintf("%.2e", mean_estimates[3]), " & ",
                        sprintf("%.2e", min_estimates[3]), " & ",
                        sprintf("%.2e", max_estimates[3]), " & ",
                        sprintf("%.2f", coverage_probabilities[3]), " & ",
                        sprintf("%.2e", rmse_estimates[3]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\alpha_{1,2}$ & ",
                        sprintf("%.2e", alpha[1, 2]), " & ",
                        sprintf("%.2e", mean_estimates[4]), " & ",
                        sprintf("%.2e", min_estimates[4]), " & ",
                        sprintf("%.2e", max_estimates[4]), " & ",
                        sprintf("%.2f", coverage_probabilities[4]), " & ",
                        sprintf("%.2e", rmse_estimates[4]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\alpha_{2,2}$ & ",
                        sprintf("%.2e", alpha[2, 2]), " & ",
                        sprintf("%.2e", mean_estimates[5]), " & ",
                        sprintf("%.2e", min_estimates[5]), " & ",
                        sprintf("%.2e", max_estimates[5]), " & ",
                        sprintf("%.2f", coverage_probabilities[5]), " & ",
                        sprintf("%.2e", rmse_estimates[5]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{1,1}$ & ",
                        sprintf("%.2e", beta[1, 1]), " & ",
                        sprintf("%.2e", mean_estimates[6]), " & ",
                        sprintf("%.2e", min_estimates[6]), " & ",
                        sprintf("%.2e", max_estimates[6]), " & ",
                        sprintf("%.2f", coverage_probabilities[6]), " & ",
                        sprintf("%.2e", rmse_estimates[6]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{2,1}$ & ",
                        sprintf("%.2e", beta[2, 1]), " & ",
                        sprintf("%.2e", mean_estimates[7]), " & ",
                        sprintf("%.2e", min_estimates[7]), " & ",
                        sprintf("%.2e", max_estimates[7]), " & ",
                        sprintf("%.2f", coverage_probabilities[7]), " & ",
                        sprintf("%.2e", rmse_estimates[7]), " \\\\", sep = "\n")



output_string <- gsub("e\\+00", "", output_string)
cat(output_string)

coverage_probabilities_1 <- colMeans(is_in_1, na.rm = TRUE)
mean_estimates_1 <- colMeans(estimates_1, na.rm = TRUE)
min_estimates_1 <- apply(estimates_1, 2, min)
max_estimates_1 <- apply(estimates_1, 2, max)
rmse_estimates_1 <- sapply(seq_len(ncol(estimates_1)), function(n) RMSE(estimates_1[, n], true_1[n]))

output_string_1 <- ""
output_string_1 <- paste0(output_string_1, "$\\beta_{1,0}$ & ",
                          sprintf("%.2e", beta0[1]), " & ",
                          sprintf("%.2e", mean_estimates_1[1]), " & ",
                          sprintf("%.2e", min_estimates_1[1]), " & ",
                          sprintf("%.2e", max_estimates_1[1]), " & ",
                          sprintf("%.2f", coverage_probabilities_1[1]), " & ",
                          sprintf("%.2e", rmse_estimates_1[1]), " \\\\", sep = "\n")

output_string_1 <- paste0(output_string_1, "$\\alpha_{1,1}$ & ",
                          sprintf("%.2e", alpha[1, 1]), " & ",
                          sprintf("%.2e", mean_estimates_1[2]), " & ",
                          sprintf("%.2e", min_estimates_1[2]), " & ",
                          sprintf("%.2e", max_estimates_1[2]), " & ",
                          sprintf("%.2f", coverage_probabilities_1[2]), " & ",
                          sprintf("%.2e", rmse_estimates_1[2]), " \\\\", sep = "\n")

output_string_1 <- paste0(output_string_1, "$\\beta_{1,1}$ & ",
                          sprintf("%.2e", beta[1, 1]), " & ",
                          sprintf("%.2e", mean_estimates_1[3]), " & ",
                          sprintf("%.2e", min_estimates_1[3]), " & ",
                          sprintf("%.2e", max_estimates_1[3]), " & ",
                          sprintf("%.2f", coverage_probabilities_1[3]), " & ",
                          sprintf("%.2e", rmse_estimates_1[3]), " \\\\", sep = "\n")


output_string_1 <- gsub("e\\+00", "", output_string_1)
cat(output_string_1)

coverage_probabilities_2 <- colMeans(is_in_2, na.rm = TRUE)
mean_estimates_2 <- colMeans(estimates_2, na.rm = TRUE)
min_estimates_2 <- apply(estimates_2, 2, min)
max_estimates_2 <- apply(estimates_2, 2, max)
rmse_estimates_2 <- sapply(seq_len(ncol(estimates_2)), function(n) RMSE(estimates_2[, n], true_2[n]))

output_string_2 <- ""
output_string_2 <- paste0(output_string_2, "$\\beta_{2,0}$ & ",
                          sprintf("%.2e", beta0[2]), " & ",
                          sprintf("%.2e", mean_estimates_2[1]), " & ",
                          sprintf("%.2e", min_estimates_2[1]), " & ",
                          sprintf("%.2e", max_estimates_2[1]), " & ",
                          sprintf("%.2f", coverage_probabilities_2[1]), " & ",
                          sprintf("%.2e", rmse_estimates_2[1]), " \\\\", sep = "\n")

output_string_2 <- paste0(output_string_2, "$\\alpha_{2,2}$ & ",
                          sprintf("%.2e", alpha[2, 2]), " & ",
                          sprintf("%.2e", mean_estimates_2[2]), " & ",
                          sprintf("%.2e", min_estimates_2[2]), " & ",
                          sprintf("%.2e", max_estimates_2[2]), " & ",
                          sprintf("%.2f", coverage_probabilities_2[2]), " & ",
                          sprintf("%.2e", rmse_estimates_2[2]), " \\\\", sep = "\n")

output_string_2 <- paste0(output_string_2, "$\\beta_{2,1}$ & ",
                          sprintf("%.2e", beta[2, 1]), " & ",
                          sprintf("%.2e", mean_estimates_2[3]), " & ",
                          sprintf("%.2e", min_estimates_2[3]), " & ",
                          sprintf("%.2e", max_estimates_2[3]), " & ",
                          sprintf("%.2f", coverage_probabilities_2[3]), " & ",
                          sprintf("%.2e", rmse_estimates_2[3]), " \\\\", sep = "\n")


output_string_2 <- gsub("e\\+00", "", output_string_2)
cat(output_string_2)
