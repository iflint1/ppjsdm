library(microbenchmark)
library(ppjsdm)
library(spatstat)

window = Rectangle_window()
window_spatstat <- owin()

configuration_spatstat <- rpoispp(lambda = 100, win = window_spatstat)
configuration <- as.Configuration(configuration_spatstat)

plot(configuration_spatstat)

set.seed(42)

b <- microbenchmark(
  "ppjsdm::gibbsm" = gibbsm(configuration, window = window, print = FALSE),
  # `logi` is the fastest method according to docs, but `mpl` appears to be faster here.
  "spatstat::ppm" = ppm(configuration_spatstat ~ 1, method = "mpl"),
  times = 1000L
)

print(b)
