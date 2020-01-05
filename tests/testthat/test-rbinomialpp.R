library(ppjsdm)
library(methods)

context("rbinomialpp")

test_that("rbinomialpp return type", {
  n <- 19
  nsim <- 14
  window <- Rectangle_window()
  sample <- rbinomialpp(window, n, nsim)

  expect_true(is.list(sample))
  expect_true(methods::is(obj = sample[[1]], class2 = "Configuration"))
  expect_equal(length(sample), nsim)
  expect_equal(length(x_coordinates(sample[[1]])), n)
})

test_that("rbinomialpp default n", {
  window <- Rectangle_window()
  set.seed(42)
  sample <- rbinomialpp(window)
  set.seed(42)
  other_sample <- rbinomialpp(window, nsim = 1)

  expect_equal(sample, other_sample)
})


test_that("rbinomialpp is in correct range (rectangle window)", {
  x_range <- c(-0.5, 4.5)
  y_range <- c(-0.8, -0.5)
  number_tests <- 100
  number_samples <- 10
  window <- Rectangle_window(x_range, y_range)

  expect_true_condition = TRUE
  for(i in seq_along(number_tests)) {
    sample <- rbinomialpp(window, n = number_samples, nsim = 1)
    x <- x_coordinates(sample)
    y <- y_coordinates(sample)
    expect_true_condition <- expect_true_condition && all(x >= x_range[1] & x <= x_range[2])
    expect_true_condition <- expect_true_condition && all(y >= y_range[1] & y <= y_range[2])
  }
  expect_true(expect_true_condition)
})

test_that("rbinomialpp is in correct range (disk window)", {
  centre <- c(-0.5, 4.5)
  radius <- 1.9
  number_tests <- 100
  number_samples <- 10
  window <- Disk_window(centre = centre, radius = radius)

  expect_true_condition = TRUE
  for(i in seq_along(number_tests)) {
    sample <- rbinomialpp(window, n = number_samples, nsim = 1)
    x <- x_coordinates(sample)
    y <- y_coordinates(sample)
    expect_true_condition <- expect_true_condition && all((x - centre[1]) * (x - centre[1]) + (y - centre[2]) * (y - centre[2]) <= radius * radius)
  }
  expect_true(expect_true_condition)
})

test_that("rbinomialpp has correct number of points.", {
  x_range <- c(-0.5, 4.5)
  y_range <- c(-0.8, -0.5)
  window <- Rectangle_window(x_range, y_range)
  n <- c(7, 2, 9, 4)

  sample <- rbinomialpp(window, n = n, nsim = 1)
  expect_equal(get_number_points(sample), c(n[1], n[2], n[3], n[4]))
  expect_equal(get_number_points(sample[1]), n[1])
  expect_equal(get_number_points(sample[2]), n[2])
  expect_equal(get_number_points(sample[3]), n[3])
  expect_equal(get_number_points(sample[4]), n[4])
})

test_that("rbinomialpp Samples from the window populated by right amount of calls to runif", {
  default_range <- c(0, 1)
  arbitrary_seed <- 42
  window <- Rectangle_window(default_range, default_range)

  set.seed(arbitrary_seed)
  expected_values <- runif(n = 4, min = 0, max = 1)
  set.seed(arbitrary_seed)
  output <- rbinomialpp(window, n = 2, nsim = 1)
  output_values <- c(x_coordinates(output), y_coordinates(output))
  expect_true(setequal(expected_values, output_values))
})
