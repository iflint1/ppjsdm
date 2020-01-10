library(microbenchmark)
library(ppjsdm)
library(spatstat)

window <-Rectangle_window()

window_spatstat <- owin()

b <- microbenchmark(
  "ppjsdm" = rppp(window = window, lambda = 10, nsim = 10),
  "spatstat" = rpoispp(10, win = window_spatstat, nsim = 10),
  times = 10000L
)

print(b)
