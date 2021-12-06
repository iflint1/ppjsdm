library(microbenchmark)
library(ppjsdm)
library(spatstat)
remove(list = ls())

set.seed(42)

b <- microbenchmark(
  "ppjsdm rstratpp" = ppjsdm::rstratpp(delta_x = 0.1, nsim = 10),
  "spatstat" = rstrat(nx = 10, nsim = 10),
  times = 1000L
)

print(b)
