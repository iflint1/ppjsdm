library(microbenchmark)
library(ppjsdm)

window <-Rectangle_window()

fast_sample <- function(n) {
  A <- matrix(c(0, 1, 0, 1), nrow = 2, ncol = 2)
  A[1,] + A[2,] * runif(2)
}

b <- microbenchmark(
  "ppjsdm" = rbinomialpp(n = 100),
  "Fast sample" = fast_sample(100),
  times = 100000L
)

print(b)
