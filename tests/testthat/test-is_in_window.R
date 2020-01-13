library(ppjsdm)

context("is_in_window")

test_that("Unique point with correct point type", {
  default_vector <- c(0, 1)
  epsilon <- 0.001
  window <- Rectangle_window(default_vector, default_vector)

  expect_true(is_in_window(0, 1, window))
  expect_true(is_in_window(1, 0, window))
  expect_true(is_in_window(0.9, 0.1, window))
  expect_true(is_in_window(1, 1, window))
  expect_false(is_in_window(-epsilon, 0.5, window))
  expect_false(is_in_window(0.5, 1 + epsilon, window))
  expect_false(is_in_window(-epsilon, 1 + epsilon, window))
})

test_that("Multiple points with correct point type", {
  default_vector <- c(0, 1)
  epsilon <- 0.001
  window <- Rectangle_window(default_vector, default_vector)

  #matrix_points_in_window <- cbind(c(0.5, 0.5), c(0.9, 0.1), c(0, 1), c(1, 0), c(1, 1))
  expect_true(is_in_window(c(0.5, 0.9, 0, 1, 1), c(0.5, 0.1, 1, 0, 1), window))

  #matrix_points_not_in_window1 <- cbind(c(0.5, 0.5), c(0.9, 0.1), c(0, 1 + epsilon), c(1, 0), c(1, 1))
  expect_false(is_in_window(c(0.5, 0.9, 0, 1, 1), c(0.5, 0.1, 1 + epsilon, 0, 1), window))

  #matrix_points_not_in_window2 <- cbind(c(0.5, 0.5), c(0.9, 0.1), c(-epsilon, 1), c(1, 0), c(1, 1))
  expect_false(is_in_window(c(0.5, 0.9, -epsilon, 1, 1), c(0.5, 0.1, 1, 0, 1), window))

  #matrix_points_not_in_window2 <- cbind(c(0.5, 0.5), c(0.9, 0.1), c(0, 1), c(1, -epsilon), c(1 + epsilon, 1))
  expect_false(is_in_window(c(0.5, 0.9, 0, 1, 1 + epsilon), c(0.5, 0.1, 1, -epsilon, 1), window))
})


test_that("Incorrect point type", {
  window <- Rectangle_window()

  expect_error(is_in_window(0, c(1, 2), window))
  expect_error(is_in_window(c(1, 2), 0, window))
})
