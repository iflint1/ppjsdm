library(microbenchmark)
library(ppjsdm)
library(spatstat)
remove(list = ls())

window_spatstat <- owin()

configuration_spatstat <- rpoispp(lambda = 100, win = window_spatstat)
configuration <- as.Configuration(configuration_spatstat)

plot(configuration_spatstat)

set.seed(42)

nd <- 100

b <- microbenchmark(
  "ppjsdm::gibbsm" = gibbsm(configuration, fitting_package = "glm", short_range = 0, dummy_factor = 1e6, max_dummy = nd * nd, print = FALSE),
  # `logi` is the fastest method according to docs, but `mpl` appears to be faster here.
  "spatstat::ppm mpl" = ppm(configuration_spatstat ~ 1, method = "mpl", nd = c(nd, nd)),
  "spatstat::ppm logi" = ppm(configuration_spatstat ~ 1, method = "logi", nd = c(nd, nd)),
  times = 100L
)

print(b)
