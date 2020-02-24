library(microbenchmark)
library(ppjsdm)
library(spatstat)

window <-Rectangle_window(c(0, 2), c(-1, 1))
window_spatstat <- owin(c(0, 2), c(-1, 1))

set.seed(42)

values <- list(c(25, 0.15), c(90, 0.08), c(15, 0.06), c(55, 0.06))

b <- lapply(values, function(v) microbenchmark(
            "ppjsdm" = rgibbs(window = window, lambda = v[1], short_range = v[2], alpha = -Inf, saturation = Inf, model = "Geyer"),
            "spatstat" = rHardcore(v[1], R = v[2], W = window_spatstat, expand = FALSE), times = 1000L))

lapply(b, function(a) print(a))
