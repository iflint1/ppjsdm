library(ppjsdm)

context("S vcov matrix computation")

test_that("compute_S works with either a sequence or a list of fits", {
  set.seed(1)
  configuration <- ppjsdm::rppp(lambda = 100)
  other <- ppjsdm::rppp(lambda = 100)
  fit <- ppjsdm::gibbsm(configuration)
  other_fit <- ppjsdm::gibbsm(other)
  expect_equal(compute_S(fit, other_fit), compute_S(list = list(fit, other_fit)))
})

test_that("compute_S works as long as fit contains theta, a regression matrix, a vector of types and a shift vector", {
  nobs <- 100
  fit_object <- NULL
  fit_object$coefficients_vector <- c(a = 1, b = 2)
  fit_object$data_list$regressors <- matrix(1, nrow = nobs, ncol = 2)
  fit_object$data_list$shift <- -5
  fit_object$data_list$type <- rep(1, nobs)

  expect_error(compute_S(fit_object), NA)
})

test_that("compute_S errors if thetas of the fits are not compatible", {
  nobs <- 100
  fit_object <- NULL
  fit_object$coefficients_vector <- c(a = 1, b = 2)
  fit_object$data_list$regressors <- matrix(1, nrow = nobs, ncol = 2)
  fit_object$data_list$shift <- -5
  fit_object$data_list$type <- rep(1, nobs)

  other_fit_object <- NULL
  other_fit_object$coefficients_vector <- c(a = 1, c = 2)
  other_fit_object$data_list$regressors <- matrix(1, nrow = nobs, ncol = 2)
  other_fit_object$data_list$shift <- -5
  other_fit_object$data_list$type <- rep(1, nobs)

  another_fit_object <- NULL
  another_fit_object$coefficients_vector <- c(a = 1, b = 2, c = 3)
  another_fit_object$data_list$regressors <- matrix(1, nrow = nobs, ncol = 2)
  another_fit_object$data_list$shift <- -5
  another_fit_object$data_list$type <- rep(1, nobs)

  expect_error(compute_S(fit_object, other_fit_object))
  expect_error(compute_S(fit_object, another_fit_object))
})
