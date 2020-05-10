library(microbenchmark)
library(ppjsdm)
library(spatstat)

window <-Rectangle_window()
window_spatstat <- owin()

set.seed(42)

values <- list(c(log(25), 0.1, 0.15), c(log(90), 0.4, 0.08), c(log(15), 0.7, 0.06), c(log(55), 1.0, 0.06))

b <- lapply(values, function(v) microbenchmark(
  "ppjsdm" = rgibbs(window = window, beta0 = v[1], short_range = v[3], alpha = log(v[2]), gamma = 0, saturation = Inf, model = "Geyer"),
  "spatstat" = rStrauss(v[1], gamma = v[2], R = v[3], W = window_spatstat, expand = FALSE), times = 1000L))

lapply(b, function(a) print(a))
