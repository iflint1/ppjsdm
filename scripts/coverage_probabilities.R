remove(list = ls())
library(ppjsdm)
library(spatstat)
set.seed(1)

window <- Rectangle_window(c(-1, 1), c(-1, 1))
nreplications <- 100
ntypes <- 2
beta0 <- c(2.5, 2)
steps <- 0
estimate_gamma <- TRUE
alpha <- cbind(c(-0.2, 0.1), c(0.1, -0.6))

covariates <- list(temperature = function(x, y) x, elevation = function(x, y) y)
beta <- cbind(c(2, 2.5), c(1, 1.5))
ndummy <- 1000

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
                          drop = FALSE,
                          saturation = saturation,
                          covariates = covariates,
                          beta = beta,
                          steps = steps)


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

library(ggplot2)
configuration <- samples[[1]]
dat <- data.frame(x = configuration$x, y = configuration$y, types = configuration$types)
X11(width = 12, height = 10)
ggplot(dat, aes(x = y, y = x)) +
  geom_point(aes(colour = types, shape = types), size = 5) +
  scale_size(breaks = seq(from = 0.16, to = 0.4, by = 0.02)) +
  coord_equal() +
  xlab("x") +
  ylab("y") +
  theme(panel.background = element_rect(fill = 'white', colour = 'red'),
        legend.title = element_blank())

X11(width = 15, height = 10)
plot(samples[[1]])

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
