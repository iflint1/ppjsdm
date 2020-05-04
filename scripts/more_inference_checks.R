library(ppjsdm)
library(spatstat)

N <- 1000
confidence <- 0.95
window <- Rectangle_window(c(0, 2), c(-1, 1))
spatstat_window <- owin(c(0, 2), c(-1, 1))

temperature <- function(x, y) x
covariates <- list(temperature = temperature)

lambda <- c(5, 5)
short_range <- cbind(c(0.1, 0.1), c(0.1, 0.1))
medium_range <- cbind(c(0.15, 0.15), c(0.15, 0.15))
long_range <- cbind(c(0.2, 0.2), c(0.2, 0.2))
alpha <- cbind(c(-1, 0.2), c(0.2, -0.5))
gamma <- cbind(c(0.8, -0.1), c(-0.1, 0.4))
beta <- cbind(c(-0.2, 0.1))
saturation <- 2
steps <- 100000
model <- "square_exponential"
medium_range_model <- "square_exponential"

configurations <- ppjsdm::rgibbs(window = window,
                                 alpha = alpha,
                                 lambda = lambda,
                                 saturation = saturation,
                                 gamma = gamma,
                                 beta = beta,
                                 covariates = covariates,
                                 model = model,
                                 medium_range_model = medium_range_model,
                                 short_range = short_range,
                                 medium_range = medium_range,
                                 long_range = long_range,
                                 nsim = N,
                                 steps = steps)
# short_range <- c(0, 0.2)
# medium_range <- c(0, 0.2)
# long_range <- c(0, 0.2)
fit <- ppjsdm::gibbsm(configurations,
                      window = window,
                      covariates = covariates,
                      model = model,
                      medium_range_model = medium_range_model,
                      short_range = short_range,
                      medium_range = medium_range,
                      long_range = long_range,
                      saturation = saturation,
                      use_glmnet = FALSE,
                      print = FALSE)
coef <- fit$coefficients
print(coef)
