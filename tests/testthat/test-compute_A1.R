library(ppjsdm)

context("Vcov matrix computation")

test_that("compute_A1 works with either a sequence or a list of fits", {
  set.seed(1)
  configuration <- ppjsdm::rppp(lambda = 100)
  other <- ppjsdm::rppp(lambda = 100)
  fit <- ppjsdm::gibbsm(configuration)
  other_fit <- ppjsdm::gibbsm(other)
  expect_equal(compute_A1(fit, other_fit), compute_A1(list = list(fit, other_fit)))
})
