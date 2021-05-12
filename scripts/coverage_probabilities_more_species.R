remove(list = ls())
library(ppjsdm)
library(spatstat)
seed <- 1

# Parameters
window <- Rectangle_window(c(0, 1), c(0, 1))
nreplications <- 10
ntypes <- 9
model <- "square_bump"
medium_range_model <- "square_exponential"
saturation <- 2
steps <- 1e5
estimate_gamma <- TRUE
diagonal_range <- seq(from = -0.6, to = 0.2, by = 0.1)
off_diagonal_range <- seq(from = -0.2, to = 0.2, by = 0.1)
beta <- cbind(seq(from = 5, to = -5, length.out = ntypes))
short_range <- matrix(0.05, ntypes, ntypes)
medium_range <- matrix(0.07, ntypes, ntypes)
long_range <- matrix(0.12, ntypes, ntypes)
covariates <- list(temperature = function(x, y) x)
max_dummy <- 2000
dummy_factor <- 1e4

# Beta0 below is the only part of the script you have to personalize if changing the number of types
# The rough idea to obtain these values is to run the script a number of times with 1 replication,
# The do sapply(seq_len(ntypes), function(n) length(samples[[1]][n]$x)) to obtain the number of individuals
# per species, and adjust beta0 so that species all have ~ same number of individuals
beta0 <- 0.9 * c(0.25, 1.75, 2.75, 5, 5, 6.5, 4.5, 5, 9)

set.seed(1)
alpha <- diag(sample(diagonal_range, size = ntypes, replace = TRUE))
alpha[lower.tri(alpha)] <- sample(off_diagonal_range, size = sum(lower.tri(alpha)), replace = TRUE)
for(i in seq_len(nrow(alpha))) {
  for(j in i:ncol(alpha)) {
    if(i != j) {
      alpha[i, j] <- alpha[j, i]
    }
  }
}

if(estimate_gamma) {
  set.seed(1)
  is_in <- matrix(NA, nrow = nreplications, ncol = ntypes + ntypes * (ntypes + 1) + length(covariates) * ntypes)
  gamma <- diag(sample(diagonal_range, size = ntypes, replace = TRUE))
  gamma[lower.tri(gamma)] <- sample(off_diagonal_range, size = sum(lower.tri(gamma)), replace = TRUE)
  for(i in seq_len(nrow(gamma))) {
    for(j in i:ncol(gamma)) {
      if(i != j) {
        gamma[i, j] <- gamma[j, i]
      }
    }
  }
} else {
  long_range <- medium_range
  is_in <- matrix(NA, nrow = nreplications, ncol = ntypes + ntypes * (ntypes + 1) / 2 + length(covariates) * ntypes)
  gamma <- matrix(0., ntypes, ntypes)
}

# The part above generated the parameters of the point process, here's the result:
print(window)
print(beta0)
print(alpha)
print(gamma)
print(short_range)
print(medium_range)
print(long_range)
print(model)
print(medium_range_model)
print(saturation)
print(covariates)
print(beta)
print(steps)

# Draw samples with those parameters
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

# Construct vector of true values
true <- beta0
for(i in seq_len(ntypes)) {
  for(j in i:ntypes) {
    true <- c(true, alpha[i, j])
  }
}
if(estimate_gamma) {
  for(i in seq_len(ntypes)) {
    for(j in i:ntypes) {
      true <- c(true, gamma[i, j])
    }
  }
}

if(length(covariates) > 0) {
  for(i in seq_len(length(covariates))) {
    for(j in seq_len(ntypes)) {
      true <- c(true, beta[j, i])
    }
  }
}

# Construct matrix in which to put estimated coefficients
estimates <- matrix(NA, nrow = nreplications, ncol = ncol(is_in))

# Fit and compute CIs
set.seed(seed)
for(i in seq_len(nreplications)) {
  fit <- gibbsm(samples[[i]],
                window = window,
                short_range = short_range,
                medium_range = medium_range,
                long_range = long_range,
                fitting_package = 'glm',
                max_dummy = max_dummy,
                dummy_factor = dummy_factor,
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

# Plot typical sample
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

# Get results in LaTeX for use in manuscript
coverage_probabilities <- colMeans(is_in, na.rm = TRUE)
mean_estimates <- colMeans(estimates, na.rm = TRUE)
median_estimates <- sapply(seq_len(ncol(estimates)), function(n) median(estimates[, n]))
mode_estimates <- sapply(seq_len(ncol(estimates)), function(n) rethinking::chainmode(estimates[, n]))
rmse_estimates <- sapply(seq_len(ncol(estimates)), function(n) sqrt(mean((estimates[, n] - true[n])^2)))

output_string <- ""
fill <- 1
for(i in seq_len(ntypes)) {
  output_string <- paste0(output_string, "$\\beta_{", i, ",0}$ & ",
                          sprintf("%.2f", beta0[i]), " & ",
                          sprintf("%.2f", mean_estimates[fill]), " & ",
                          sprintf("%.2f", rmse_estimates[fill]), " & ",
                          sprintf("%.2f", coverage_probabilities[fill]), " \\\\", sep = "\n")
  fill <- fill + 1
}

for(i in seq_len(ntypes)) {
  for(j in i:ntypes) {
    output_string <- paste0(output_string, "$\\alpha_{", i, ",", j, "}$ & ",
                            sprintf("%.2f", alpha[i, j]), " & ",
                            sprintf("%.2f", mean_estimates[fill]), " & ",
                            sprintf("%.2f", rmse_estimates[fill]), " & ",
                            sprintf("%.2f", coverage_probabilities[fill]), " \\\\", sep = "\n")
    fill <- fill + 1
  }
}

if(estimate_gamma) {
  for(i in seq_len(ntypes)) {
    for(j in i:ntypes) {
      output_string <- paste0(output_string, "$\\gamma_{", i, ",", j, "}$ & ",
                              sprintf("%.2f", gamma[i, j]), " & ",
                              sprintf("%.2f", mean_estimates[fill]), " & ",
                              sprintf("%.2f", rmse_estimates[fill]), " & ",
                              sprintf("%.2f", coverage_probabilities[fill]), " \\\\", sep = "\n")
      fill <- fill + 1
    }
  }
}

if(length(covariates) > 0) {
  for(i in seq_len(length(covariates))) {
    for(j in seq_len(ntypes)) {
      output_string <- paste0(output_string, "$\\beta_{", j, ",", i, "}$ & ",
                              sprintf("%.2f", beta[j, i]), " & ",
                              sprintf("%.2f", mean_estimates[fill]), " & ",
                              sprintf("%.2f", rmse_estimates[fill]), " & ",
                              sprintf("%.2f", coverage_probabilities[fill]), " \\\\", sep = "\n")
      fill <- fill + 1
    }
  }
}

output_string <- gsub("e\\+00", "", output_string)
cat(output_string)

# Average estimates
average_output_string <- ""
fill <- 0
average_output_string <- paste0(average_output_string, "$\\beta_{0}$ & ",
                                sprintf("%.2f", mean(true[seq_len(ntypes)])), " & ",
                                sprintf("%.2f", mean(mean_estimates[seq_len(ntypes)])), " & ",
                                sprintf("%.2f", mean(rmse_estimates[seq_len(ntypes)])), " & ",
                                sprintf("%.2f", mean(coverage_probabilities[seq_len(ntypes)])), " \\\\", sep = "\n")
fill <- fill + ntypes

average_output_string <- paste0(average_output_string, "$\\alpha$ & ",
                                sprintf("%.2f", mean(true[fill + seq_len(ntypes * (ntypes + 1) / 2)])), " & ",
                                sprintf("%.2f", mean(mean_estimates[fill + seq_len(ntypes * (ntypes + 1) / 2)])), " & ",
                                sprintf("%.2f", mean(rmse_estimates[fill + seq_len(ntypes * (ntypes + 1) / 2)])), " & ",
                                sprintf("%.2f", mean(coverage_probabilities[fill + seq_len(ntypes * (ntypes + 1) / 2)])), " \\\\", sep = "\n")

fill <- fill + ntypes * (ntypes + 1) / 2

if(estimate_gamma) {
  average_output_string <- paste0(average_output_string, "$\\gamma$ & ",
                                  sprintf("%.2f", mean(true[fill + seq_len(ntypes * (ntypes + 1) / 2)])), " & ",
                                  sprintf("%.2f", mean(mean_estimates[fill + seq_len(ntypes * (ntypes + 1) / 2)])), " & ",
                                  sprintf("%.2f", mean(rmse_estimates[fill + seq_len(ntypes * (ntypes + 1) / 2)])), " & ",
                                  sprintf("%.2f", mean(coverage_probabilities[fill + seq_len(ntypes * (ntypes + 1) / 2)])), " \\\\", sep = "\n")
  fill <- fill + ntypes * (ntypes + 1) / 2
}

if(length(covariates) > 0) {
  for(i in seq_len(length(covariates))) {
    average_output_string <- paste0(average_output_string, "$\\beta_{", i, "}$ & ",
                                    sprintf("%.2f", mean(true[fill + seq_len(ntypes)])), " & ",
                                    sprintf("%.2f", mean(mean_estimates[fill + seq_len(ntypes)])), " & ",
                                    sprintf("%.2f", mean(rmse_estimates[fill + seq_len(ntypes)])), " & ",
                                    sprintf("%.2f", mean(coverage_probabilities[fill + seq_len(ntypes)])), sep = "\n")
    fill <- fill + ntypes
  }
}

average_output_string <- gsub("e\\+00", "", average_output_string)
cat(average_output_string)
