library(ppjsdm)

context("Configuration")

test_that("Correct result in simple case", {
  configuration <- Configuration(c(0.5, 1.5, 2.5), c(0.5, 0.5, 0.5), factor(c("a", "b", "b")))

  expect_equal(compute_delta_phi_dispersion(configuration, c(1, 0.5), 0, 2, "i", 0), c(-0.5, 0.25))
  expect_equal(compute_delta_phi_dispersion(configuration, c(1, 0.5), 2, 3, "i", 0), c(-0.5, -1.0, 0.0))
})
