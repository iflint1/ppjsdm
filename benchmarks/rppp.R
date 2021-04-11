library(microbenchmark)
library(ppjsdm)
library(spatstat)
remove(list = ls())

window <- Rectangle_window()
window_spatstat <- owin()

set.seed(42)

b <- microbenchmark(
  "ppjsdm rppp" = rppp(window = window, lambda = 10, nsim = 10),
  "ppjsdm rgibbs" = rgibbs(window = window, beta0 = log(10), nsim = 10),
  "spatstat" = rpoispp(10, win = window_spatstat, nsim = 10),
  times = 10000L
)

print(b)
