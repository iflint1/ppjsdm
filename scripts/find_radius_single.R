remove(list = ls())
library(ggplot2)
library(ppjsdm)
library(spatstat)
seed <- 1

window <- Rectangle_window(c(0, 1), c(0, 1))
nreplications <- 1e3
ntypes <- 2
beta0 <- c(6.5, 2.6)
steps <- 1e6
alpha <- cbind(c(-1, -0.5), c(-0.5, 2))

covariates <- list(x = function(x, y) x - 0.5)
beta <- cbind(c(-1, 1))
max_dummy <- 1e4
dummy_factor <- 2

model <- "exponential"
medium_range_model <- "square_exponential"
saturation <- 2

short_range_search <- seq(from = 0.0025, to = 0.5, by = 0.0025)

short_range <- cbind(c(0.05, 0.05), c(0.05, 0.05))
medium_range <- short_range
long_range <- medium_range
gamma <- matrix(0., ntypes, ntypes)

set.seed(seed)
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
dat <- data.frame(x = configuration$x,
                  y = configuration$y,
                  types = configuration$types)
png(file = "unknown_radius.png", bg = "white", width = 600, height = 400)
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

tm <- Sys.time()

results1 <- list()
for(k1 in seq_len(length(short_range_search))) {
  test_short_range <- matrix(short_range_search[k1], ncol = 2, nrow = 2)
  average_logLik <- mean(sapply(seq_len(nreplications), function(i) {
    set.seed(seed)
    fit <- ppjsdm::gibbsm(samples[[i]],
                          window = window,
                          short_range = test_short_range,
                          medium_range = medium_range,
                          long_range = long_range,
                          fitting_package = "glm",
                          max_dummy = max_dummy,
                          dummy_factor = dummy_factor,
                          model = model,
                          medium_range_model = medium_range_model,
                          covariates = covariates,
                          saturation = saturation)
    logLik(fit$complete)
  }))
  results1[[length(results1) + 1]] <- list(short_range = test_short_range, average_logLik = average_logLik)
}

Sys.time() - tm

tm <- Sys.time()

results2 <- lapply(seq_len(nreplications), function(i) {
  lapply(seq_len(length(short_range_search)), function(k1) {
    test_short_range <- matrix(short_range_search[k1], ncol = 2, nrow = 2)
    set.seed(seed)
    fit <- ppjsdm::gibbsm(samples[[i]],
                          window = window,
                          short_range = test_short_range,
                          medium_range = medium_range,
                          long_range = long_range,
                          fitting_package = "glm",
                          max_dummy = max_dummy,
                          dummy_factor = dummy_factor,
                          model = model,
                          medium_range_model = medium_range_model,
                          covariates = covariates,
                          saturation = saturation)
    list(logLik = logLik(fit$complete),
         short = short_range_search[k1])
  })
})

Sys.time() - tm

best_short <- lapply(results1, function(res) res$short_range)[[which.max(sapply(results1, function(res) res$average_logLik))]]
print(best_short)

df <- data.frame(short_range = short_range_search, average_logLik = sapply(results1, function(res) res$average_logLik))
ggplot(data = df, aes(x = short_range, y = average_logLik)) +
  geom_line(size = 0.5) +
  geom_point() +
  geom_vline(xintercept = short_range[1, 1], colour = 'red') +
  xlab("Short-range interaction radius") +
  ylab("Average pseudo-loglikelihood") +
  theme_minimal()

best_shorts <- sapply(results2, function(r) r[[which.max(sapply(r, function(rr) rr$logLik))]]$short)

print(mean(best_shorts[best_shorts < max(short_range_search)]))
print(sum(best_shorts < max(short_range_search)) / length(best_shorts))

df <- data.frame(x = best_shorts)
ggplot(data = df, aes(x = x)) +
  geom_histogram(binwidth = 0.005) +
  geom_vline(xintercept = short_range[1, 1], colour = 'red') +
  geom_vline(xintercept = median(best_shorts), colour = 'blue') +
  geom_vline(xintercept = mean(best_shorts), colour = 'yellow') +
  xlab("Short-range interaction radius") +
  theme_minimal()

which(best_shorts == max(short_range_search))

plot(x = short_range_search, y = sapply(results2[[88]], function(a) a$logLik), type = 'l')


library(ggplot2)
configuration <- samples[[1]]
dat <- data.frame(x = configuration$x, y = configuration$y, types = configuration$types)
X11(width = 12, height = 10)
ggplot(dat, aes(x = x, y = y)) +
  geom_point(aes(colour = types, shape = types), size = 5) +
  scale_size(breaks = seq(from = 0.16, to = 0.4, by = 0.02)) +
  coord_equal() +
  xlab("x") +
  ylab("y") +
  theme(panel.background = element_rect(fill = 'white', colour = 'red'),
        legend.title = element_blank())

tm <- Sys.time()


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
