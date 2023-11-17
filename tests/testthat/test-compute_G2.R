library(ppjsdm)

context("G2 vcov matrix computation")

test_that("compute_G2 works with either a sequence or a list of fits", {
  set.seed(1)
  configuration <- ppjsdm::rppp(lambda = 100)
  other <- ppjsdm::rppp(lambda = 100)
  fit <- ppjsdm::gibbsm(configuration)
  other_fit <- ppjsdm::gibbsm(other)

  set.seed(1)
  result1 <- ppjsdm::compute_G2(fit, other_fit)

  set.seed(1)
  result2 <- ppjsdm::compute_G2(list = list(fit, other_fit))
  expect_equal(result1, result2)
})
