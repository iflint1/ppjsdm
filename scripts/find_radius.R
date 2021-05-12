remove(list = ls())
library(ppjsdm)
library(spatstat)
seed <- 1

window <- Rectangle_window(c(0, 1), c(0, 1))
nreplications <- 10
ntypes <- 2
beta0 <- c(4, 5)
steps <- 0
alpha <- cbind(c(-2, -1), c(-1, -2))

covariates <- list(x = function(x, y) x - 0.5)
beta <- cbind(c(-1, 1))
max_dummy <- 1000
dummy_factor <- 4

model <- "exponential"
medium_range_model <- "square_exponential"
saturation <- 2

short_range_search <- seq(from = 0.02, to = 0.3, by = 0.01)

short_range <- cbind(c(0.02, 0.08), c(0.08, 0.2))
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

results <- list()
for(k1 in seq_len(length(short_range_search))) {
  for(k2 in seq_len(length(short_range_search))) {
    for(k3 in seq_len(length(short_range_search))) {
      test_short_range <- cbind(c(short_range_search[k1], short_range_search[k2]),
                                c(short_range_search[k2], short_range_search[k3]))
      average_logLik <- mean(sapply(seq_len(nreplications), function(i) {
        set.seed(seed)
        fit <- ppjsdm::gibbsm(samples[[i]],
                              window = window,
                              short_range = test_short_range,
                              medium_range = medium_range,
                              long_range = long_range,
                              fitting_package = 'glm',
                              max_dummy = max_dummy,
                              dummy_factor = dummy_factor,
                              model = model,
                              medium_range_model = medium_range_model,
                              covariates = covariates,
                              saturation = saturation)
        logLik(fit$complete)
      }))
      results[[length(results) + 1]] <- list(short_range = test_short_range, average_logLik = average_logLik)
    }
  }
}

# Code below optimizes radii for each replication

# best_short <- lapply(seq_len(nreplications), function(i) {
#   best_logLik <- -Inf
#   best_short <- NA
#   for(k1 in seq_len(length(short_range_search))) {
#     for(k2 in seq_len(length(short_range_search))) {
#       for(k3 in seq_len(length(short_range_search))) {
#         test_short_range <- cbind(c(short_range_search[k1], short_range_search[k2]), c(short_range_search[k2], short_range_search[k3]))
#         set.seed(seed)
#         fit <- ppjsdm::gibbsm(samples[[i]],
#                               window = window,
#                               short_range = test_short_range,
#                               medium_range = medium_range,
#                               long_range = long_range,
#                               fitting_package = 'glm',
#                               max_dummy = max_dummy,
#                               dummy_factor = dummy_factor,
#                               model = model,
#                               medium_range_model = medium_range_model,
#                               covariates = covariates,
#                               saturation = saturation)
#         current_logLik <- logLik(fit$complete)
#         if(current_logLik > best_logLik) {
#           best_logLik <- current_logLik
#           best_short <- test_short_range
#         }
#       }
#     }
#   }
#
#   best_short
# })

Sys.time() - tm

plot(sapply(results, function(res) res$average_logLik))

best_short <- lapply(results, function(res) res$short_range)[[which.max(sapply(results, function(res) res$average_logLik))]]

# print(Reduce("+", best_short) / length(best_short))
print(best_short)

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
