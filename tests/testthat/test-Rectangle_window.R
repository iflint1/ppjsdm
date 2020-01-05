library(ppjsdm)

context("Rectangle_window")

test_that("Rectangle_window is an S4 class", {
  window <- Rectangle_window()

  expect_is(window, "Rectangle_window")
})

test_that("Rectangle_window accessors", {
  x_coordinates <- c(0, 1)
  y_coordinates <- c(2, 3)
  window <- Rectangle_window(x_coordinates, y_coordinates)

  expect_true(is.vector(x_range(window)))
  expect_true(all.equal(x_range(window), x_coordinates))

  expect_true(is.vector(y_range(window)))
  expect_true(all.equal(y_range(window), y_coordinates))
})


test_that("Rectangle_window default constructor", {
  default_range <- c(0, 1)
  default_window <- Rectangle_window()

  expect_identical(x_range(default_window), default_range)
  expect_identical(y_range(default_window), default_range)
})

test_that("Rectangle_window specified arguments constructor", {
  vector1 <- c(-1, 2)
  vector2 <- c(0.5, 1.5)
  window_specified_arguments <- Rectangle_window(y_range = vector2, x_range = vector1)

  expect_identical(x_range(window_specified_arguments), vector1)
  expect_identical(y_range(window_specified_arguments), vector2)
})

test_that("Rectangle_window unspecified arguments constructors", {
  vector1 <- c(-1, 2)
  vector2 <- c(0.5, 1.5)
  window_unspecified_arguments <- Rectangle_window(vector1, vector2)

  expect_identical(x_range(window_unspecified_arguments), vector1)
  expect_identical(y_range(window_unspecified_arguments), vector2)
})

test_that("Rectangle_window one argument constructor", {
  default_range <- c(0, 1)
  vector1 <- c(-1, 2)
  window_one_argument <- Rectangle_window(vector1)

  expect_identical(x_range(window_one_argument), vector1)
  expect_identical(y_range(window_one_argument), default_range)
})

test_that("Rectangle_window copy constructor", {
  vector1 <- c(-1, 2)
  vector2 <- c(0.5, 1.5)
  window_specified_arguments <- Rectangle_window(x_range = vector1, y_range = vector2)
  copied_window <- Rectangle_window(window_specified_arguments)

  expect_identical(x_range(window_specified_arguments), x_range(copied_window))
  expect_identical(y_range(window_specified_arguments), y_range(copied_window))
})

test_that("Rectangle_window unsupported constructors", {
  default_range <- c(0, 1)

  expect_error(Rectangle_window(c(0, 1, 3)))
  expect_error(Rectangle_window(c('a', 'b')))
  expect_error(Rectangle_window(c(1, 0)))
  expect_error(Rectangle_window(default_range, c(0, 1, 3)))
  expect_error(Rectangle_window(default_range, c('a', 'b')))
  expect_error(Rectangle_window(default_range, c(1, 0)))
})
