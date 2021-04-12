library(microbenchmark)
library(ppjsdm)
library(spatstat)
remove(list = ls())

window_spatstat <- owin()

temperature <- function(x, y) x + y

configuration_spatstat <- rpoispp(lambda = function(x, y) exp(4 + temperature(x, y)), win = window_spatstat)
configuration <- as.Configuration(configuration_spatstat)

plot(configuration_spatstat)

set.seed(42)

nd <- 100

b <- microbenchmark(
  "ppjsdm::gibbsm" = gibbsm(configuration, covariates = list(temperature = temperature), use_glmnet = FALSE, short_range = 0, dummy_factor = 1e6, max_dummy = nd * nd, print = FALSE),
  # `logi` is the fastest method according to docs, but `mpl` appears to be faster here.
  "spatstat::ppm" = ppm(configuration_spatstat ~ 1 + temperature, method = "mpl", nd = c(nd, nd)),
  times = 100L
)

print(b)
