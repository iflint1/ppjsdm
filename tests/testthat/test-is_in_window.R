library(ppjsdm)

context("is_in_window")

test_that("is_in_window with unique point with correct point type", {
  default_vector <- c(0, 1)
  epsilon <- 0.001
  window <- Rectangle_window(default_vector, default_vector)

  expect_true(is_in_window(c(0, 1), window))
  expect_true(is_in_window(c(1, 0), window))
  expect_true(is_in_window(c(0.9, 0.1), window))
  expect_true(is_in_window(c(1, 1), window))
  expect_false(is_in_window(c(-epsilon, 0.5), window))
  expect_false(is_in_window(c(0.5, 1 + epsilon), window))
  expect_false(is_in_window(c(-epsilon, 1 + epsilon), window))
})

test_that("is_in_window with multiple points with correct point type", {
  default_vector <- c(0, 1)
  epsilon <- 0.001
  window <- Rectangle_window(default_vector, default_vector)

  matrix_points_in_window <- cbind(c(0.5, 0.5), c(0.9, 0.1), c(0, 1), c(1, 0), c(1, 1))
  expect_true(is_in_window(matrix_points_in_window, window))

  matrix_points_not_in_window1 <- cbind(c(0.5, 0.5), c(0.9, 0.1), c(0, 1 + epsilon), c(1, 0), c(1, 1))
  expect_false(is_in_window(matrix_points_not_in_window1, window))

  matrix_points_not_in_window2 <- cbind(c(0.5, 0.5), c(0.9, 0.1), c(-epsilon, 1), c(1, 0), c(1, 1))
  expect_false(is_in_window(matrix_points_not_in_window2, window))

  matrix_points_not_in_window2 <- cbind(c(0.5, 0.5), c(0.9, 0.1), c(0, 1), c(1, -epsilon), c(1 + epsilon, 1))
  expect_false(is_in_window(matrix_points_not_in_window2, window))
})


test_that("is_in_window with incorrect point type", {
  window <- Rectangle_window()

  #expect_error(is_in_window(c('a', 'b'), window))
  expect_error(is_in_window(c(0), window))
  expect_error(is_in_window(0.5, window))
  expect_error(is_in_window(c(0, 1, 2), window))
})
