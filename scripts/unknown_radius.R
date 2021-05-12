remove(list = ls())
library(ppjsdm)
library(spatstat)
set.seed(1)

window <- Rectangle_window(c(0, 1), c(0, 1))
nreplications <- 100
ntypes <- 2
beta0 <- c(4, 3.5)
steps <- 100000
alpha <- cbind(c(0.4, -0.3), c(-0.3, 0.4))

covariates <- list(temperature = function(x, y) x)
beta <- cbind(c(1.5, 2))
ndummy <- 1000

model <- "square_bump"
medium_range_model <- "square_exponential"
saturation <- 2

short_range <- cbind(c(0.04, 0.06), c(0.06, 0.03))
wrong_short_range1 <- short_range - 0.02
wrong_short_range2 <- short_range + 0.02
medium_range <- cbind(c(0.07, 0.02), c(0.02, 0.06))
long_range <- medium_range
gamma <- matrix(0., ntypes, ntypes)
is_in1 <- matrix(NA, nrow = nreplications, ncol = ntypes + ntypes * (ntypes + 1) / 2 + length(covariates) * ntypes)
estimates1 <- matrix(NA, nrow = nreplications, ncol = ncol(is_in1))
is_in2 <- matrix(NA, nrow = nreplications, ncol = ncol(is_in1))
estimates2 <- matrix(NA, nrow = nreplications, ncol = ncol(estimates1))
true <- c(beta0, alpha[1, 1], alpha[1, 2], alpha[2, 2])

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

tm <- Sys.time()

for(i in seq_len(nreplications)) {
  fit1 <- gibbsm(samples[[i]],
                 window = window,
                 short_range = wrong_short_range1,
                 medium_range = medium_range,
                 long_range = long_range,
                 fitting_package = 'glm',
                 ndummy = ndummy,
                 model = model,
                 medium_range_model = medium_range_model,
                 covariates = covariates,
                 saturation = saturation)

  fit2 <- gibbsm(samples[[i]],
                 window = window,
                 short_range = wrong_short_range2,
                 medium_range = medium_range,
                 long_range = long_range,
                 fitting_package = 'glm',
                 ndummy = ndummy,
                 model = model,
                 medium_range_model = medium_range_model,
                 covariates = covariates,
                 saturation = saturation)

  s1 <- summary(fit1)$coefficients
  estimate1 <- s1$coefficients

  lower1 <- s1$CI95_lo
  upper1 <- s1$CI95_hi
  is_in1[i, ] <- true >= lower1 & true <= upper1
  estimates1[i, ] <- estimate1

  s2 <- summary(fit2)$coefficients
  estimate2 <- s2$coefficients

  lower2 <- s2$CI95_lo
  upper2 <- s2$CI95_hi
  is_in2[i, ] <- true >= lower2 & true <= upper2
  estimates2[i, ] <- estimate2
}

Sys.time() - tm

coverage_probabilities1 <- colMeans(is_in1, na.rm = TRUE)
mean_estimates1 <- colMeans(estimates1, na.rm = TRUE)
coverage_probabilities2 <- colMeans(is_in2, na.rm = TRUE)
mean_estimates2 <- colMeans(estimates2, na.rm = TRUE)

output_string <- ""
output_string <- paste0(output_string, "$\\beta_{1,0}$ & ",
                        sprintf("%.2e", beta0[1]), " & ",
                        sprintf("%.2e", mean_estimates1[1]), " & ",
                        sprintf("%.2f", coverage_probabilities1[1]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{2,0}$ & ",
                        sprintf("%.2e", beta0[2]), " & ",
                        sprintf("%.2e", mean_estimates1[2]), " & ",
                        sprintf("%.2f", coverage_probabilities1[2]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\alpha_{1,1}$ & ",
                        sprintf("%.2e", alpha[1, 1]), " & ",
                        sprintf("%.2e", mean_estimates1[3]), " & ",
                        sprintf("%.2f", coverage_probabilities1[3]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\alpha_{1,2}$ & ",
                        sprintf("%.2e", alpha[1, 2]), " & ",
                        sprintf("%.2e", mean_estimates1[4]), " & ",
                        sprintf("%.2f", coverage_probabilities1[4]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\alpha_{2,2}$ & ",
                        sprintf("%.2e", alpha[2, 2]), " & ",
                        sprintf("%.2e", mean_estimates1[5]), " & ",
                        sprintf("%.2f", coverage_probabilities1[5]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{1,1}$ & ",
                        sprintf("%.2e", beta[1, 1]), " & ",
                        sprintf("%.2e", mean_estimates1[6]), " & ",
                        sprintf("%.2f", coverage_probabilities1[6]), " \\\\", sep = "\n")

output_string <- paste0(output_string, "$\\beta_{2,1}$ & ",
                        sprintf("%.2e", beta[2, 1]), " & ",
                        sprintf("%.2e", mean_estimates1[7]), " & ",
                        sprintf("%.2f", coverage_probabilities1[7]), " \\\\", sep = "\n")



output_string <- gsub("e\\+00", "", output_string)
cat(output_string)

output_string2 <- ""
output_string2 <- paste0(output_string2, "$\\beta_{1,0}$ & ",
                        sprintf("%.2e", beta0[1]), " & ",
                        sprintf("%.2e", mean_estimates2[1]), " & ",
                        sprintf("%.2f", coverage_probabilities2[1]), " \\\\", sep = "\n")

output_string2 <- paste0(output_string2, "$\\beta_{2,0}$ & ",
                        sprintf("%.2e", beta0[2]), " & ",
                        sprintf("%.2e", mean_estimates2[2]), " & ",
                        sprintf("%.2f", coverage_probabilities2[2]), " \\\\", sep = "\n")

output_string2 <- paste0(output_string2, "$\\alpha_{1,1}$ & ",
                        sprintf("%.2e", alpha[1, 1]), " & ",
                        sprintf("%.2e", mean_estimates2[3]), " & ",
                        sprintf("%.2f", coverage_probabilities2[3]), " \\\\", sep = "\n")

output_string2 <- paste0(output_string2, "$\\alpha_{1,2}$ & ",
                        sprintf("%.2e", alpha[1, 2]), " & ",
                        sprintf("%.2e", mean_estimates2[4]), " & ",
                        sprintf("%.2f", coverage_probabilities2[4]), " \\\\", sep = "\n")

output_string2 <- paste0(output_string2, "$\\alpha_{2,2}$ & ",
                        sprintf("%.2e", alpha[2, 2]), " & ",
                        sprintf("%.2e", mean_estimates2[5]), " & ",
                        sprintf("%.2f", coverage_probabilities2[5]), " \\\\", sep = "\n")

output_string2 <- paste0(output_string2, "$\\beta_{1,1}$ & ",
                        sprintf("%.2e", beta[1, 1]), " & ",
                        sprintf("%.2e", mean_estimates2[6]), " & ",
                        sprintf("%.2f", coverage_probabilities2[6]), " \\\\", sep = "\n")

output_string2 <- paste0(output_string2, "$\\beta_{2,1}$ & ",
                        sprintf("%.2e", beta[2, 1]), " & ",
                        sprintf("%.2e", mean_estimates2[7]), " & ",
                        sprintf("%.2f", coverage_probabilities2[7]), " \\\\", sep = "\n")



output_string2 <- gsub("e\\+00", "", output_string2)
cat(output_string2)
