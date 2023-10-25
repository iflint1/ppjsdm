library(ppjsdm)

context("Rectangle_window")

test_that("Rectangle_window has correct type", {
  expect_is(ppjsdm::Rectangle_window(), "Rectangle_window")
})

test_that("Rectangle_window range", {
  x_coordinates <- cumsum(stats::runif(n = 2))
  y_coordinates <- cumsum(stats::runif(n = 2))
  window <- ppjsdm::Rectangle_window(x_range = x_coordinates, y_range = y_coordinates)

  expect_true(is.vector(ppjsdm::x_range(window)))
  expect_true(all.equal(ppjsdm::x_range(window), x_coordinates))

  expect_true(is.vector(ppjsdm::y_range(window)))
  expect_true(all.equal(ppjsdm::y_range(window), y_coordinates))
})

test_that("Rectangle_window default constructor", {
  default_range <- c(0, 1)
  default_window <- ppjsdm::Rectangle_window()

  expect_identical(default_window, ppjsdm::Rectangle_window(default_range, default_range))
})

test_that("Rectangle_window inverted interval", {
  range <- c(-1, 2)
  inverted_range <- c(2, -1)
  window <- ppjsdm::Rectangle_window(range)
  inverted_window <- ppjsdm::Rectangle_window(inverted_range)

  expect_identical(window, inverted_window)
})

test_that("Rectangle_window one argument constructor", {
  default_range <- c(0, 1)
  vector1 <- cumsum(stats::runif(n = 2))
  window_one_argument <- ppjsdm::Rectangle_window(vector1)

  expect_identical(window_one_argument, ppjsdm::Rectangle_window(x_range = vector1, y_range = default_range))
})

test_that("Rectangle_window copy constructor", {
  vector1 <- cumsum(stats::runif(n = 2))
  vector2 <- cumsum(stats::runif(n = 2))
  window_specified_arguments <- ppjsdm::Rectangle_window(x_range = vector1, y_range = vector2)
  copied_window <- ppjsdm::Rectangle_window(window_specified_arguments)

  expect_identical(window_specified_arguments, copied_window)
})

test_that("Rectangle_window unsupported constructors", {
  default_range <- c(0, 1)

  expect_error(ppjsdm::Rectangle_window(c(0, 1, 3)))
  expect_error(ppjsdm::Rectangle_window(c('a', 'b')))
  expect_error(ppjsdm::Rectangle_window(default_range, c(0, 1, 3)))
  expect_error(ppjsdm::Rectangle_window(default_range, c('a', 'b')))
})
