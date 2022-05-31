library(ppjsdm)
library(spatstat)

set.seed(1)

alpha <- cbind(c(0.3, -0.2), c(-0.2, -0.4))
gamma <- cbind(c(-0.3, 0.2), c(0.2, 0.4))
beta0 <- c(4.5, 4)
short_range <- matrix(0.05, 2, 2)
medium_range <- matrix(0.1, 2, 2)
long_range <- matrix(0.2, 2, 2)
steps <- 1e5
type <- 2
nsamples <- 100

under_fitted_alpha <- alpha
under_fitted_alpha[2, 2] <- alpha[2, 2] - 0.5
under_fitted_beta0 <- beta0 + c(0, 0.4) # Obtained by trail and error

over_fitted_alpha <- alpha
over_fitted_alpha[2, 2] <- alpha[2, 2] + 0.5
over_fitted_beta0 <- beta0 - c(0, 0.6) # Obtained by trail and error

under_fitted_gamma <- gamma
under_fitted_gamma[2, 2] <- gamma[2, 2] - 0.5
medium_under_fitted_beta0 <- beta0 + c(0, 0.7) # Obtained by trail and error

over_fitted_gamma <- gamma
over_fitted_gamma[2, 2] <- gamma[2, 2] + 0.5
medium_over_fitted_beta0 <- beta0 - c(0, 0.9) # Obtained by trail and error

configuration <- ppjsdm::rgibbs(alpha = alpha, beta0 = beta0, gamma = gamma, short_range = short_range,
                                medium_range = medium_range, long_range = long_range, steps = steps)
plot(configuration)

#
# set.seed(1)
# configuration <- ppjsdm::rgibbs(alpha = alpha, beta0 = medium_over_fitted_beta0, gamma = over_fitted_gamma, short_range = short_range,
#                                 medium_range = medium_range, long_range = long_range, steps = steps)
# configuration[1]
# configuration[2]
# plot(configuration)


set.seed(1)
conditional_samples <- lapply(seq_len(nsamples), function(i) {
  z1 <- Configuration(x = configuration$x[configuration$types %in% levels(configuration$types)[-type]],
                      y = configuration$y[configuration$types %in% levels(configuration$types)[-type]],
                      types = configuration$types[configuration$types %in% levels(configuration$types)[-type]])
  z2 <- Configuration(x = configuration$x[configuration$types %in% levels(configuration$types)[type]],
                      y = configuration$y[configuration$types %in% levels(configuration$types)[type]],
                      types = configuration$types[configuration$types %in% levels(configuration$types)[type]])
  result <- ppjsdm::rgibbs(alpha = alpha,
                           beta0 = beta0,
                           gamma = gamma,
                           short_range = short_range,
                           medium_range = medium_range,
                           long_range = long_range,
                           # starting_configuration = z2,
                           steps = steps / 10,
                           only_simulate_these_types = type,
                           conditional_configuration = z1)
})

set.seed(1)
under_conditional_samples <- lapply(seq_len(nsamples), function(i) {
    z1 <- Configuration(x = configuration$x[configuration$types %in% levels(configuration$types)[-type]],
                        y = configuration$y[configuration$types %in% levels(configuration$types)[-type]],
                        types = configuration$types[configuration$types %in% levels(configuration$types)[-type]])
    z2 <- Configuration(x = configuration$x[configuration$types %in% levels(configuration$types)[type]],
                        y = configuration$y[configuration$types %in% levels(configuration$types)[type]],
                        types = configuration$types[configuration$types %in% levels(configuration$types)[type]])
    ppjsdm::rgibbs(alpha = under_fitted_alpha,
                   beta0 = under_fitted_beta0,
                   gamma = gamma,
                   short_range = short_range,
                   medium_range = medium_range,
                   long_range = long_range,
                   # starting_configuration = z2,
                   steps = steps / 10,
                   only_simulate_these_types = type,
                   conditional_configuration = z1)
})

set.seed(1)
over_conditional_samples <- lapply(seq_len(nsamples), function(i) {
  z1 <- Configuration(x = configuration$x[configuration$types %in% levels(configuration$types)[-type]],
                      y = configuration$y[configuration$types %in% levels(configuration$types)[-type]],
                      types = configuration$types[configuration$types %in% levels(configuration$types)[-type]])
  z2 <- Configuration(x = configuration$x[configuration$types %in% levels(configuration$types)[type]],
                      y = configuration$y[configuration$types %in% levels(configuration$types)[type]],
                      types = configuration$types[configuration$types %in% levels(configuration$types)[type]])
  ppjsdm::rgibbs(alpha = over_fitted_alpha,
                 beta0 = over_fitted_beta0,
                 gamma = gamma,
                 short_range = short_range,
                 medium_range = medium_range,
                 long_range = long_range,
                 # starting_configuration = z2,
                 steps = steps / 10,
                 only_simulate_these_types = type,
                 conditional_configuration = z1)
})

set.seed(1)
medium_under_conditional_samples <- lapply(seq_len(nsamples), function(i) {
  z1 <- Configuration(x = configuration$x[configuration$types %in% levels(configuration$types)[-type]],
                      y = configuration$y[configuration$types %in% levels(configuration$types)[-type]],
                      types = configuration$types[configuration$types %in% levels(configuration$types)[-type]])
  z2 <- Configuration(x = configuration$x[configuration$types %in% levels(configuration$types)[type]],
                      y = configuration$y[configuration$types %in% levels(configuration$types)[type]],
                      types = configuration$types[configuration$types %in% levels(configuration$types)[type]])
  ppjsdm::rgibbs(alpha = alpha,
                 beta0 = medium_under_fitted_beta0,
                 gamma = under_fitted_gamma,
                 short_range = short_range,
                 medium_range = medium_range,
                 long_range = long_range,
                 # starting_configuration = z2,
                 steps = steps / 10,
                 only_simulate_these_types = type,
                 conditional_configuration = z1)
})

set.seed(1)
medium_over_conditional_samples <- lapply(seq_len(nsamples), function(i) {
  z1 <- Configuration(x = configuration$x[configuration$types %in% levels(configuration$types)[-type]],
                      y = configuration$y[configuration$types %in% levels(configuration$types)[-type]],
                      types = configuration$types[configuration$types %in% levels(configuration$types)[-type]])
  z2 <- Configuration(x = configuration$x[configuration$types %in% levels(configuration$types)[type]],
                      y = configuration$y[configuration$types %in% levels(configuration$types)[type]],
                      types = configuration$types[configuration$types %in% levels(configuration$types)[type]])
  ppjsdm::rgibbs(alpha = alpha,
                 beta0 = medium_over_fitted_beta0,
                 gamma = over_fitted_gamma,
                 short_range = short_range,
                 medium_range = medium_range,
                 long_range = long_range,
                 # starting_configuration = z2,
                 steps = steps / 10,
                 only_simulate_these_types = type,
                 conditional_configuration = z1)
})

plot(conditional_samples[[1]])
plot(configuration)

xseq <- seq(from = 0, to = 1, length.out = 100)
grid <- expand.grid(x = xseq, y = xseq)

result <- ppjsdm::compute_papangelou(configuration = configuration[-type],
                                     x = grid[, 1],
                                     y = grid[, 2],
                                     type = type,
                                     mark = 1.0,
                                     alpha = alpha,
                                     beta0 = beta0,
                                     gamma = gamma,
                                     short_range = short_range,
                                     medium_range = medium_range,
                                     long_range = long_range)
papangelou <- log(spatstat.geom::im(matrix(result,
                                           nrow = length(xseq),
                                           ncol = length(xseq),
                                           byrow = TRUE),
                                    xrange = c(0, 1),
                                    yrange = c(0, 1)))
plot(papangelou)
points(configuration[type])

ggplot_envelope <- function(env, horizontal_line_at_one = TRUE, ylim) {
  require(ggplot2)

  df <- data.frame(r = env$r,
                   obs = env$obs,
                   mean = env$mmean,
                   lo = env$lo,
                   hi = env$hi)[-1, ]

  if(!missing(ylim)) {
    df$hi <- pmin(df$hi, ylim[2])
  }

  g <- ggplot(df) +
    theme_minimal() +
    geom_line(aes(x = r, y = obs), colour = "black") +
    geom_line(aes(x = r, y = mean), colour = "red", linetype = "dashed") +
    geom_ribbon(aes(x = r, ymin = lo, ymax = hi), fill = "grey", alpha = 0.5) +
    xlab("r") +
    ylab("")
  if(!missing(ylim)) {
    g <- g + ylim(ylim)
  }

  if(horizontal_line_at_one) {
    g <- g + geom_line(aes(x = r, y = rep(1, nrow(df))), colour = "green", alpha = 0.5)
  }
  g
}

ggplot_envelope(ylim = c(0, 10), envelope(fun = pcfinhom,
                                          # lambda = (papangelou),
                                          Y = as.ppp(configuration[type], W = owin()),
                                          normpower = 2,
                                          nsim = length(over_conditional_samples),
                                          r = seq(from = 0, to = 0.25, length.out = 1e3),
                                          simulate = lapply(over_conditional_samples, function(conf) as.ppp(conf, W = owin()))))

ggplot_envelope(ylim = c(0, 10), envelope(fun = pcfinhom,
                                          # lambda = (papangelou),
                                          Y = as.ppp(configuration[type], W = owin()),
                                          normpower = 2,
                                          nsim = length(conditional_samples),
                                          r = seq(from = 0, to = 0.25, length.out = 1e3),
                                          simulate = lapply(conditional_samples, function(conf) as.ppp(conf[type], W = owin()))))

ggplot_envelope(ylim = c(0, 10), envelope(fun = pcfinhom,
                                          # lambda = (papangelou),
                                          Y = as.ppp(configuration[type], W = owin()),
                                          normpower = 2,
                                          nsim = length(under_conditional_samples),
                                          r = seq(from = 0, to = 0.25, length.out = 1e3),
                                          simulate = lapply(under_conditional_samples, function(conf) as.ppp(conf[type], W = owin()))))

ggplot_envelope(ylim = c(0, 10), envelope(fun = pcfinhom,
                                          # lambda = (papangelou),
                                          Y = as.ppp(configuration[type], W = owin()),
                                          normpower = 2,
                                          nsim = length(medium_under_conditional_samples),
                                          r = seq(from = 0, to = 0.25, length.out = 1e3),
                                          simulate = lapply(medium_under_conditional_samples, function(conf) as.ppp(conf[type], W = owin()))))

ggplot_envelope(ylim = c(0, 10), envelope(fun = pcfinhom,
                                          # lambda = (papangelou),
                                          Y = as.ppp(configuration[type], W = owin()),
                                          normpower = 2,
                                          nsim = length(medium_over_conditional_samples),
                                          r = seq(from = 0, to = 0.25, length.out = 1e3),
                                          simulate = lapply(medium_over_conditional_samples, function(conf) as.ppp(conf[type], W = owin()))))




ggplot_envelope(ylim = c(0, 10), envelope(fun = pcf,
                                          Y = as.ppp(configuration[type], W = owin()),
                                          nsim = length(over_conditional_samples),
                                          r = seq(from = 0, to = 0.25, length.out = 1e3),
                                          simulate = lapply(over_conditional_samples, function(conf) as.ppp(conf, W = owin()))))

ggplot_envelope(ylim = c(0, 10), envelope(fun = pcf,
                                          Y = as.ppp(configuration[type], W = owin()),
                                          nsim = length(conditional_samples),
                                          r = seq(from = 0, to = 0.25, length.out = 1e3),
                                          simulate = lapply(conditional_samples, function(conf) as.ppp(conf[type], W = owin()))))

ggplot_envelope(ylim = c(0, 10), envelope(fun = pcf,
                                          Y = as.ppp(configuration[type], W = owin()),
                                          nsim = length(under_conditional_samples),
                                          r = seq(from = 0, to = 0.25, length.out = 1e3),
                                          simulate = lapply(under_conditional_samples, function(conf) as.ppp(conf[type], W = owin()))))

ggplot_envelope(ylim = c(0, 10), envelope(fun = pcf,
                                          Y = as.ppp(configuration[type], W = owin()),
                                          nsim = length(medium_under_conditional_samples),
                                          r = seq(from = 0, to = 0.25, length.out = 1e3),
                                          simulate = lapply(medium_under_conditional_samples, function(conf) as.ppp(conf[type], W = owin()))))

ggplot_envelope(ylim = c(0, 10), envelope(fun = pcf,
                                          Y = as.ppp(configuration[type], W = owin()),
                                          nsim = length(medium_over_conditional_samples),
                                          r = seq(from = 0, to = 0.25, length.out = 1e3),
                                          simulate = lapply(medium_over_conditional_samples, function(conf) as.ppp(conf[type], W = owin()))))

