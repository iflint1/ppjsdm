library(microbenchmark)
library(ppjsdm)
library(spatstat)

window <-Rectangle_window()
window_spatstat <- owin()

set.seed(42)

b <- microbenchmark(
  "ppjsdm" = rgibbs(window = window, lambda = 25, alpha = -Inf, saturation = Inf, model = "Geyer", short_range = 0.1),
  "spatstat" = rHardcore(25, R = 0.1, W = window_spatstat, expand = FALSE),
  times = 10000L
)

print(b)
