library(microbenchmark)
library(ppjsdm)
library(spatstat)

window <-Rectangle_window()
window_spatstat <- owin()

set.seed(42)

values <- list(c(log(25), 0.15), c(log(90), 0.08), c(log(15), 0.06), c(log(55), 0.06))

b <- lapply(values, function(v) microbenchmark(
            "ppjsdm" = rgibbs(window = window, beta0 = v[1], short_range = v[2], alpha = -Inf, gamma = 0, saturation = Inf, model = "Geyer"),
            "spatstat" = rHardcore(v[1], R = v[2], W = window_spatstat, expand = FALSE), times = 1000L))

lapply(b, function(a) print(a))
