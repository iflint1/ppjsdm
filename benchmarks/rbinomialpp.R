library(microbenchmark)
library(ppjsdm)

window <-Rectangle_window()

fast_sample <- function(n) {
  A <- matrix(c(0, 1, 0, 1), nrow = 2, ncol = 2)
  A[1,] + A[2,] * runif(2)
  #runif(2 * n)
}

microbenchmark(
  "Default" = rbinomialpp(window, 10),
  "Fast sample" = fast_sample(10),
  times = 100000L
)
