library(microbenchmark)
library(ppjsdm)
library(spatstat)
remove(list = ls())

window <- Rectangle_window()
window_spatstat <- owin()

set.seed(42)

nthreads <- 4
values <- list(c(log(25), 0.1, 0.15), c(log(90), 0.4, 0.08), c(log(15), 0.7, 0.06), c(log(55), 1.0, 0.06))

# CFTP
lapply(values, function(v) microbenchmark(
  "ppjsdm" = rgibbs(window = window, nthreads = nthreads, beta0 = v[1], alpha = 0.5 * log(v[2]),
                    short_range = v[3], gamma = 0, saturation = Inf, model = "Geyer", nsim = 10),
  "spatstat" = rStrauss(exp(v[1]), gamma = v[2], R = v[3], nsim = 10, W = window_spatstat, expand = FALSE), times = 100L))

# Metropolis-Hastings
nsteps <- 1e4
lapply(values, function(v) microbenchmark(
  "ppjsdm" = rgibbs(window = window, nthreads = nthreads, beta0 = v[1], alpha = 0.5 * log(v[2]),
                    short_range = v[3], gamma = 0, saturation = Inf, model = "Geyer", steps = nsteps, nsim = 10),
  "spatstat" = {
    mod01 <- list(cif = "strauss", par = list(beta = exp(v[1]), gamma = v[2], r = v[3]),
                  w = window_spatstat)
    rmh(model = mod01, nsim = 10, control = list(nrep = nsteps), verbose = FALSE)
  }, times = 10L))
