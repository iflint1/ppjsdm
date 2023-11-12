# The goal of this script is to provide a series of test cases for ppjsdm potting functions.
# This script should be executed every time the functions are modified to make sure
# all the plots look as expected.

library(ppjsdm)

set.seed(1)

ntypes <- 4
alpha <- matrix(0, ntypes, ntypes)
covariates <- list(x = function(x, y) x - 0.5, y = function(x, y) y - 0.5)
short_range <- list(matrix(0.05, ntypes, ntypes), matrix(0.1, ntypes, ntypes))
medium_range <- matrix(0.15, ntypes, ntypes)
long_range <- matrix(0.2, ntypes, ntypes)

# The three samples below have the same average number of points, but the types labels differ
samples1 <- ppjsdm::rppp(lambda = setNames(rep(10, ntypes), nm = paste0("t", seq_len(ntypes))))
samples2 <- ppjsdm::rppp(lambda = setNames(rep(10, ntypes), nm = paste0("t", 1 + seq_len(ntypes))))
samples3 <- ppjsdm::rppp(lambda = setNames(rep(10, ntypes), nm = paste0("t", 2 + seq_len(ntypes))))

# Fits are constructed below
fit1 <- ppjsdm::gibbsm(samples1,
                       covariates = covariates,
                       short_range = short_range,
                       medium_range = medium_range,
                       long_range = long_range)
fit2 <- ppjsdm::gibbsm(samples2,
                       covariates = covariates[1], # This fit uses only one covariate
                       short_range = short_range,
                       medium_range = medium_range,
                       long_range = long_range)
fit3 <- ppjsdm::gibbsm(samples3,
                       covariates = covariates,
                       short_range = short_range[1], # This fit uses only one potential
                       medium_range = medium_range,
                       long_range = long_range)

# First, run a few tests to make sure the summaries are computed automatically
plot(fit1, fit2, fit3)
box_plot(fit1, fit2, fit3, coefficient = "alpha2")
plot(fit1, list = list(fit2, fit3), coefficient = "beta")
box_plot(A = fit1, list = list(C = fit2, fit3), coefficient = "y")

chord_diagram_plot(fit1, only_statistically_significant = FALSE)
chord_diagram_plot(fit1, only_statistically_significant = FALSE, coefficient = "gamma")

# Below, we rewrite over the fitted coefficients and the corresponding summaries to make
# it easier to inspect the plots. All the fitted coefficients have a straightforward dependence on
# the species indices.
# Compute the summaries
summ <- list(summary(fit1), summary(fit2), summary(fit3))

# In order to help check the plots, we standardise the fitted coefficients
for(i in seq_len(ntypes)) {
  fit1$coefficients$beta[i, 1] <- i
  fit2$coefficients$beta[i, 1] <- 2 * i
  fit3$coefficients$beta[i, 1] <- -i

  fit1$coefficients$beta[i, 2] <- 2 * i
  fit3$coefficients$beta[i, 2] <- 0.5 * i

  summ[[1]]$lo$beta[i, 1] <- i - 1
  summ[[1]]$hi$beta[i, 1] <- i + 1
  summ[[1]]$lo_numerical$beta[i, 1] <- i - 1 / 2
  summ[[1]]$hi_numerical$beta[i, 1] <- i + 1 / 2

  summ[[2]]$lo$beta[i, 1] <- 2 * (i - 1)
  summ[[2]]$hi$beta[i, 1] <- 2 * (i + 1)
  summ[[2]]$lo_numerical$beta[i, 1] <- 2 * (i - 1 / 2)
  summ[[2]]$hi_numerical$beta[i, 1] <- 2 * (i + 1 / 2)

  summ[[3]]$lo$beta[i, 1] <- -i - 1
  summ[[3]]$hi$beta[i, 1] <- -i + 1
  summ[[3]]$lo_numerical$beta[i, 1] <- -i - 1 / 2
  summ[[3]]$hi_numerical$beta[i, 1] <- -i + 1 / 2

  summ[[1]]$lo$beta[i, 2] <- 2 * (i - 1)
  summ[[1]]$hi$beta[i, 2] <- 2 * (i + 1)
  summ[[1]]$lo_numerical$beta[i, 2] <- 2 * (i - 1 / 2)
  summ[[1]]$hi_numerical$beta[i, 2] <- 2 * (i + 1 / 2)

  summ[[3]]$lo$beta[i, 2] <- 0.5 * (i - 1)
  summ[[3]]$hi$beta[i, 2] <- 0.5 * (i + 1)
  summ[[3]]$lo_numerical$beta[i, 2] <- 0.5 * (i - 1 / 2)
  summ[[3]]$hi_numerical$beta[i, 2] <- 0.5 * (i + 1 / 2)

  for(j in seq_len(ntypes)) {
    fit1$coefficients$gamma[i, j] <- i + j
    fit2$coefficients$gamma[i, j] <- 0.5 * (i + j)
    fit3$coefficients$gamma[i, j] <- (i + j - 3 * abs(i - j))

    fit1$coefficients$alpha[[1]][i, j] <- abs(i - j)
    fit2$coefficients$alpha[[1]][i, j] <- 0.5 * abs(i - j)
    fit3$coefficients$alpha[[1]][i, j] <- -abs(i - j)

    fit1$coefficients$alpha[[2]][i, j] <- abs(i - 2 * j)
    fit2$coefficients$alpha[[2]][i, j] <- 0.5 * abs(i - 2 * j)

    summ[[1]]$lo$gamma[i, j] <- i + j - 1
    summ[[1]]$hi$gamma[i, j] <- i + j + 1
    summ[[1]]$lo_numerical$gamma[i, j] <- i + j - 1 / 2
    summ[[1]]$hi_numerical$gamma[i, j] <- i + j + 1 / 2

    summ[[2]]$lo$gamma[i, j] <- 0.5 * (i + j - 1)
    summ[[2]]$hi$gamma[i, j] <- 0.5 * (i + j + 1)
    summ[[2]]$lo_numerical$gamma[i, j] <- 0.5 * (i + j - 1 / 2)
    summ[[2]]$hi_numerical$gamma[i, j] <- 0.5 * (i + j + 1 / 2)

    summ[[3]]$lo$gamma[i, j] <- (i + j - 3 * abs(i - j)) - 1
    summ[[3]]$hi$gamma[i, j] <- (i + j - 3 * abs(i - j)) + 1
    summ[[3]]$lo_numerical$gamma[i, j] <- (i + j - 3 * abs(i - j)) - 1 / 2
    summ[[3]]$hi_numerical$gamma[i, j] <- (i + j - 3 * abs(i - j)) + 1 / 2

    summ[[1]]$lo$alpha[[1]][i, j] <- abs(i - j) - 1
    summ[[1]]$hi$alpha[[1]][i, j] <- abs(i - j) + 1
    summ[[1]]$lo_numerical$alpha[[1]][i, j] <- abs(i - j) - 1 / 2
    summ[[1]]$hi_numerical$alpha[[1]][i, j] <- abs(i - j) + 1 / 2

    summ[[2]]$lo$alpha[[1]][i, j] <- 0.5 * (abs(i - j) - 1)
    summ[[2]]$hi$alpha[[1]][i, j] <- 0.5 * (abs(i - j) + 1)
    summ[[2]]$lo_numerical$alpha[[1]][i, j] <- 0.5 * (abs(i - j) - 1 / 2)
    summ[[2]]$hi_numerical$alpha[[1]][i, j] <- 0.5 * (abs(i - j) + 1 / 2)

    summ[[3]]$lo$alpha[[1]][i, j] <- -abs(i - j) - 1
    summ[[3]]$hi$alpha[[1]][i, j] <- -abs(i - j) + 1
    summ[[3]]$lo_numerical$alpha[[1]][i, j] <- -abs(i - j) - 1 / 2
    summ[[3]]$hi_numerical$alpha[[1]][i, j] <- -abs(i - j) + 1 / 2

    summ[[1]]$lo$alpha[[2]][i, j] <- abs(i - 2 * j) - 1
    summ[[1]]$hi$alpha[[2]][i, j] <- abs(i - 2 * j) + 1
    summ[[1]]$lo_numerical$alpha[[2]][i, j] <- abs(i - 2 * j) - 1 / 2
    summ[[1]]$hi_numerical$alpha[[2]][i, j] <- abs(i - 2 * j) + 1 / 2

    summ[[2]]$lo$alpha[[2]][i, j] <- 0.5 * (abs(i - 2 * j) - 1)
    summ[[2]]$hi$alpha[[2]][i, j] <- 0.5 * (abs(i - 2 * j) + 1)
    summ[[2]]$lo_numerical$alpha[[2]][i, j] <- 0.5 * (abs(i - 2 * j) - 1 / 2)
    summ[[2]]$hi_numerical$alpha[[2]][i, j] <- 0.5 * (abs(i - 2 * j) + 1 / 2)
  }
}

print(fit1$coefficients$alpha[[1]])
#    t1 t2 t3 t4
# t1  0  1  2  3
# t2  1  0  1  2
# t3  2  1  0  1
# t4  3  2  1  0
print(fit2$coefficients$alpha[[1]])
#     t2  t3  t4  t5
# t2 0.0 0.5 1.0 1.5
# t3 0.5 0.0 0.5 1.0
# t4 1.0 0.5 0.0 0.5
# t5 1.5 1.0 0.5 0.0
print(fit3$coefficients$alpha[[1]])
#    t3 t4 t5 t6
# t3  0 -1 -2 -3
# t4 -1  0 -1 -2
# t5 -2 -1  0 -1
# t6 -3 -2 -1  0

print(fit1$coefficients$alpha[[2]])
#    t1 t2 t3 t4
# t1  1  3  5  7
# t2  0  2  4  6
# t3  1  1  3  5
# t4  2  0  2  4
print(fit2$coefficients$alpha[[2]])
#     t2  t3  t4  t5
# t2 0.5 1.5 2.5 3.5
# t3 0.0 1.0 2.0 3.0
# t4 0.5 0.5 1.5 2.5
# t5 1.0 0.0 1.0 2.0

print(fit1$coefficients$gamma)
#    t1 t2 t3 t4
# t1  2  3  4  5
# t2  3  4  5  6
# t3  4  5  6  7
# t4  5  6  7  8
print(fit2$coefficients$gamma)
#     t2  t3  t4  t5
# t2 1.0 1.5 2.0 2.5
# t3 1.5 2.0 2.5 3.0
# t4 2.0 2.5 3.0 3.5
# t5 2.5 3.0 3.5 4.0
print(fit3$coefficients$gamma)
#    t3 t4 t5 t6
# t3  2  0 -2 -4
# t4  0  4  2  0
# t5 -2  2  6  4
# t6 -4  0  4  8

print(fit1$coefficients$beta)
#    x y
# t1 1 2
# t2 2 4
# t3 3 6
# t4 4 8
print(fit2$coefficients$beta)
#    x
# t2 2
# t3 4
# t4 6
# t5 8
print(fit3$coefficients$beta)
#     x   y
# t3 -1 0.5
# t4 -2 1.0
# t5 -3 1.5
# t6 -4 2.0

plot(fit1, fit2, fit3, summ = summ, coefficient = "alpha")
plot(fit1, fit2, fit3, summ = summ, coefficient = "alpha1")
plot(fit1, fit2, fit3, summ = summ, coefficient = "alpha2")
plot(fit1, fit2, fit3, summ = summ, coefficient = "gamma")
plot(fit1, fit2, fit3, summ = summ, coefficient = "beta")
plot(fit1, fit2, fit3, summ = summ, coefficient = c("y", "x"))
plot(fit1, fit2, fit3, summ = summ, coefficient = "beta1")
plot(fit1, fit2, fit3, summ = summ, coefficient = "beta2")
plot(fit1, fit2, fit3, summ = summ, coefficient = "x")
plot(fit1, fit2, fit3, summ = summ, coefficient = "y")

plot(fit3, summ = summ[3], coefficient = "gamma", classes = c(t1 = "A", t2 = "A", type3 = "B", t4 = "B"),
     full_names = c(t3 = "type3", t4 = "type4"))

chord_diagram_plot(fit1, summ = summ[[1]], only_statistically_significant = FALSE, coefficient = "alpha1")

chord_diagram_plot(fit3, summ = summ[[3]], only_statistically_significant = FALSE, coefficient = "alpha1")

chord_diagram_plot(fit3, summ = summ[[3]], compute_confidence_intervals = FALSE,
                   coefficient = "alpha")

chord_diagram_plot(fit1, summ = summ[[1]], only_statistically_significant = FALSE,
                   coefficient = "alpha", classes = c(t1 = "A", t2 = "A", t3 = "B", t4 = "B"))

chord_diagram_plot(fit1, summ = summ[[1]], only_statistically_significant = FALSE,
                   coefficient = "alpha", classes = c(t1 = "A", t2 = "B", t3 = "C", t4 = "D", t5 = "E", t5 = "E"))

chord_diagram_plot(fit1, summ = summ[[1]], include_self = FALSE,
                   only_statistically_significant = FALSE,
                   coefficient = "alpha", classes = c(t1 = "A", t2 = "A", t3 = "B", t4 = "B"))

chord_diagram_plot(fit1, summ = summ[[1]], show_grid_ticks = FALSE,
                   only_statistically_significant = FALSE,
                   coefficient = "alpha", classes = c(t1 = "A", t2 = "A", t3 = "B", t4 = "B"))

chord_diagram_plot(fit1, summ = summ[[1]], show_grid_ticks = FALSE,
                   only_statistically_significant = FALSE,
                   coefficient = "alpha1", classes = c(t1 = "A", t2 = "A", t3 = "B", t4 = "B"))

chord_diagram_plot(fit3, summ = summ[[3]], only_statistically_significant = FALSE,
                   coefficient = "gamma", classes = c(t1 = "A", t2 = "A", t3 = "B", t4 = "B"))

chord_diagram_plot(fit3, summ = summ[[3]], only_statistically_significant = FALSE,
                   coefficient = "gamma")

chord_diagram_plot(fit3, summ = summ[[3]], ninteractions = 5, only_statistically_significant = FALSE,
                   coefficient = "gamma")

chord_diagram_plot(fit3, summ = summ[[3]], classes = c(t1 = "A", t2 = "A", t3 = "B", t4 = "B"),
                   ninteractions = 2, only_statistically_significant = FALSE,
                   coefficient = "gamma")

chord_diagram_plot(fit3, summ = summ[[3]], only_statistically_significant = FALSE,
                   coefficient = "gamma", classes = c(t1 = "A", t2 = "A", t3 = "B", t4 = "B"),
                   full_names = c(t3 = "type3", t4 = "type4"))

chord_diagram_plot(fit3, summ = summ[[3]], only_statistically_significant = FALSE,
                   coefficient = "gamma", classes = c(t1 = "A", t2 = "A", type3 = "B", t4 = "B"),
                   full_names = c(t3 = "type3", t4 = "type4"))

chord_diagram_plot(fit3, summ = summ[[3]], only_statistically_significant = FALSE,
                   coefficient = "gamma", classes = c(t1 = "A", t2 = "A", type3 = "B", t4 = "B"),
                   full_names = c(t3 = "type3", t4 = "type4"), outward_facing_names = TRUE)

chord_diagram_plot(fit3, summ = summ[[3]], show_grid_ticks = FALSE,
                   only_statistically_significant = FALSE,
                   coefficient = "gamma", classes = c(t1 = "A", t2 = "A", type3 = "B", t4 = "B"),
                   full_names = c(t3 = "type3", t4 = "type4"), outward_facing_names = TRUE)

chord_diagram_plot(fit3, summ = summ[[3]], show_grid_ticks = FALSE,
                   only_statistically_significant = FALSE,
                   coefficient = "gamma", classes = c(t1 = "A", t2 = "A", type3 = "B", t4 = "B"),
                   full_names = c(t3 = "type3", t4 = "type4"), outward_facing_names = TRUE, circle_margin = 1)

