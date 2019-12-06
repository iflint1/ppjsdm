library(ppjsdm)

context("rppp")

test_that("rppp return type", {
  n <- 5
  lambda <- 1.5
  window <- Rectangle_window()
  sample <- rppp(window, lambda = lambda, nsim = n)

  expect_true(is.list(sample))
  expect_true(methods::is(obj = sample[[1]], class2 = "Configuration"))
  expect_equal(length(sample), n)
})

setClass(Class = "Bad_window", slots = list(alpha = "numeric"))

test_that("rppp wrong window type", {
  lambda <- 1.5
  window <- methods::new("Bad_window", alpha = 0)

  expect_error(rppp(window, lambda))
})

test_that("rppp default n", {
  lambda <- 1.5
  window <- Rectangle_window()
  sample <- rppp(window, lambda)

  expect_equal(length(sample), 1)
})

test_that("rppp default lambda", {
  n <- 6
  window <- Rectangle_window()

  set.seed(42)
  sample <- rppp(window, nsim = 6)
  set.seed(42)
  other_sample <- rppp(window, lambda = 1, nsim = 6)

  expect_true(all.equal(sample, other_sample))
})


test_that("rppp points are in correct range", {
  x_range <- c(-0.5, 4.5)
  y_range <- c(-0.8, -0.5)
  lambda <- 1.5
  number_tests <- 100
  window <- Rectangle_window(x_range, y_range)

  expect_true_condition = TRUE
  for(i in seq_along(number_tests)) {
    sample <- rppp(window, lambda, nsim = 1)
    x <- x_coordinates(sample)
    y <- y_coordinates(sample)
    expect_true_condition <- expect_true_condition && all(x >= x_range[1] & x <= x_range[2])
    expect_true_condition <- expect_true_condition && all(y >= y_range[1] & y <= y_range[2])
  }
  expect_true(expect_true_condition)
})

test_that("rppp Number of points populated by correct number of calls to rpois", {
  lambda <- 150
  arbitrary_seed <- 42
  default_range <- c(0, 1)
  window <- Rectangle_window(default_range, default_range)

  set.seed(arbitrary_seed)
  expected_values <- rpois(1, lambda)
  set.seed(arbitrary_seed)
  output <- rppp(window, lambda, nsim = 1)
  output_values <- length(x_coordinates(output))
  expect_equal(expected_values, output_values)
})
