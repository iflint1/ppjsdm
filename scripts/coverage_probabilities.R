remove(list = ls())
library(ppjsdm)
library(spatstat)
seed <- 1
set.seed(seed)

window <- Rectangle_window(c(-1, 1), c(-1, 1))
nreplications <- 1e3
ntypes <- 2
beta0 <- c(2.5, 2)
steps <- 0
estimate_gamma <- TRUE
alpha <- cbind(c(-0.2, 0.1), c(0.1, -0.6))

covariates <- list(temperature = function(x, y) x, elevation = function(x, y) y)
beta <- cbind(c(2, 2.5), c(1, 1.5))
max_dummy <- 2000
dummy_factor <- 1e5

model <- "square_bump"
medium_range_model <- "square_exponential"
saturation <- 2

short_range <- matrix(0.05, ntypes, ntypes)
medium_range <- matrix(0.07, ntypes, ntypes)
if(estimate_gamma) {
  long_range <- matrix(0.12, ntypes, ntypes)
  is_in <- matrix(NA, nrow = nreplications, ncol = ntypes + ntypes * (ntypes + 1) + length(covariates) * ntypes)
  gamma <- cbind(c(-0.6, -0.3), c(-0.3, 0.))
} else {
  long_range <- matrix(0.07, ntypes, ntypes)
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
                          saturation = saturation,
                          covariates = covariates,
                          beta = beta,
                          drop = FALSE,
                          steps = steps)

library(ggplot2)
configuration <- samples[[1]]
dat <- data.frame(x = configuration$x,
                  y = configuration$y,
                  types = configuration$types)
png(file = "first_experiment_typical_sample.png", width = 600, height = 400)
ggplot(dat, aes(x = x, y = y)) +
  geom_point(aes(colour = types, shape = types), size = 4) +
  scale_size(breaks = seq(from = 0.16, to = 0.4, by = 0.02)) +
  coord_equal() +
  theme_minimal(base_size = 30) +
  xlab("") +
  ylab("") +
  xlim(x_range(window)) +
  ylim(y_range(window)) +
  guides(color = guide_legend(override.aes = list(size = 8))) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 40),
        axis.text = element_text(colour = "black"))
dev.off()

set.seed(seed)
for(i in seq_len(nreplications)) {
  fit <- gibbsm(samples[[i]],
                window = window,
                short_range = short_range,
                medium_range = medium_range,
                long_range = long_range,
                max_dummy = max_dummy,
                dummy_factor = dummy_factor,
                use_regularization = FALSE,
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

coverage_probabilities <- colMeans(is_in, na.rm = TRUE)
mean_estimates <- colMeans(estimates, na.rm = TRUE)
median_estimates <- sapply(seq_len(ncol(estimates)), function(n) median(estimates[, n]))
mode_estimates <- sapply(seq_len(ncol(estimates)), function(n) rethinking::chainmode(estimates[, n]))
rmse_estimates <- sapply(seq_len(ncol(estimates)), function(n) sqrt(mean((estimates[, n] - true[n])^2)))

output_string <- ""
output_string <- paste0(output_string, "$\\beta_{1,0}$ & ",
                        sprintf("%.2e", beta0[1]), " & ",
                        sprintf("%.2e", mean_estimates[1]), " & ",
                        sprintf("%.2e", median_estimates[1]), " & ",
                        sprintf("%.2e", rmse_estimates[1]), " & ",
                        sprintf("%.2f", coverage_probabilities[1]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{2,0}$ & ",
                        sprintf("%.2e", beta0[2]), " & ",
                        sprintf("%.2e", mean_estimates[2]), " & ",
                        sprintf("%.2e", median_estimates[2]), " & ",
                        sprintf("%.2e", rmse_estimates[2]), " & ",
                        sprintf("%.2f", coverage_probabilities[2]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\alpha_{1,1}$ & ",
                        sprintf("%.2e", alpha[1, 1]), " & ",
                        sprintf("%.2e", mean_estimates[3]), " & ",
                        sprintf("%.2e", median_estimates[3]), " & ",
                        sprintf("%.2e", rmse_estimates[3]), " & ",
                        sprintf("%.2f", coverage_probabilities[3]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\alpha_{1,2}$ & ",
                        sprintf("%.2e", alpha[1, 2]), " & ",
                        sprintf("%.2e", mean_estimates[4]), " & ",
                        sprintf("%.2e", median_estimates[4]), " & ",
                        sprintf("%.2e", rmse_estimates[4]), " & ",
                        sprintf("%.2f", coverage_probabilities[4]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\alpha_{2,2}$ & ",
                        sprintf("%.2e", alpha[2, 2]), " & ",
                        sprintf("%.2e", mean_estimates[5]), " & ",
                        sprintf("%.2e", median_estimates[5]), " & ",
                        sprintf("%.2e", rmse_estimates[5]), " & ",
                        sprintf("%.2f", coverage_probabilities[5]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\gamma_{1,1}$ & ",
                        sprintf("%.2e", gamma[1, 1]), " & ",
                        sprintf("%.2e", mean_estimates[6]), " & ",
                        sprintf("%.2e", median_estimates[6]), " & ",
                        sprintf("%.2e", rmse_estimates[6]), " & ",
                        sprintf("%.2f", coverage_probabilities[6]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\gamma_{1,2}$ & ",
                        sprintf("%.2e", gamma[1, 2]), " & ",
                        sprintf("%.2e", mean_estimates[7]), " & ",
                        sprintf("%.2e", median_estimates[7]), " & ",
                        sprintf("%.2e", rmse_estimates[7]), " & ",
                        sprintf("%.2f", coverage_probabilities[7]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\gamma_{2,2}$ & ",
                        sprintf("%.2e", gamma[2, 2]), " & ",
                        sprintf("%.2e", mean_estimates[8]), " & ",
                        sprintf("%.2e", median_estimates[8]), " & ",
                        sprintf("%.2e", rmse_estimates[8]), " & ",
                        sprintf("%.2f", coverage_probabilities[8]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{1,1}$ & ",
                        sprintf("%.2e", beta[1, 1]), " & ",
                        sprintf("%.2e", mean_estimates[9]), " & ",
                        sprintf("%.2e", median_estimates[9]), " & ",
                        sprintf("%.2e", rmse_estimates[9]), " & ",
                        sprintf("%.2f", coverage_probabilities[9]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{2,1}$ & ",
                        sprintf("%.2e", beta[2, 1]), " & ",
                        sprintf("%.2e", mean_estimates[10]), " & ",
                        sprintf("%.2e", median_estimates[10]), " & ",
                        sprintf("%.2e", rmse_estimates[10]), " & ",
                        sprintf("%.2f", coverage_probabilities[10]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{1,2}$ & ",
                        sprintf("%.2e", beta[1, 2]), " & ",
                        sprintf("%.2e", mean_estimates[11]), " & ",
                        sprintf("%.2e", median_estimates[11]), " & ",
                        sprintf("%.2e", rmse_estimates[11]), " & ",
                        sprintf("%.2f", coverage_probabilities[11]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{2,2}$ & ",
                        sprintf("%.2e", beta[2, 2]), " & ",
                        sprintf("%.2e", mean_estimates[12]), " & ",
                        sprintf("%.2e", median_estimates[12]), " & ",
                        sprintf("%.2e", rmse_estimates[12]), " & ",
                        sprintf("%.2f", coverage_probabilities[12]), " \\\\", sep = "\n")



output_string <- gsub("e\\+00", "", output_string)
cat(output_string)
