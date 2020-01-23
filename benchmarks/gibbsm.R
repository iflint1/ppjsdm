library(microbenchmark)
library(ppjsdm)
library(spatstat)

window <- Rectangle_window()
window_spatstat <- owin()

configuration <- rppp(lambda = 100, window = window)
configuration_spatstat <- rpoispp(lambda = 100, win = window_spatstat)

set.seed(42)

b <- microbenchmark(
  "ppjsdm::gibbsm" = gibbsm(configuration, window = window, print = FALSE),
  # `logi` is the fastest method according to docs, but `mpl` appears to be faster here.
  "spatstat::ppm" = ppm(configuration_spatstat ~ 1, method = "mpl"),
  times = 1000L
)

print(b)
