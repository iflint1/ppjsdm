library(microbenchmark)
library(ppjsdm)
library(spatstat)
remove(list = ls())

window <- Rectangle_window()
window_spatstat <- owin()

set.seed(42)

values <- list(c(log(25), 0.15), c(log(90), 0.08), c(log(15), 0.06), c(log(55), 0.06))

lapply(values, function(v) microbenchmark(
  "ppjsdm" = rgibbs(window = window, beta0 = v[1], short_range = v[2], alpha = -Inf, gamma = 0, saturation = Inf, model = "Geyer", nsim = 10),
  "spatstat" = rHardcore(exp(v[1]), R = v[2], W = window_spatstat, expand = FALSE, nsim = 10), times = 100L))
