library(ppjsdm)

context("rppp")

test_that("Return type.", {
  n <- 5
  sample <- rppp(nsim = n)

  expect_true(is.list(sample))
  expect_true(methods::is(obj = sample[[1]], class2 = "Configuration"))
  expect_equal(length(sample), n)

  sample <- rppp(nsim = 1, drop = TRUE)

  expect_true(methods::is(obj = sample, class2 = "Configuration"))
})

test_that("Wrong window type.", {
  window <- structure(list(x_range = x_range, y_range = y_range))
  expect_error(rppp(window = window))
})

test_that("Default arguments.", {
  set.seed(42)
  sample <- rppp()
  set.seed(42)
  other_sample <- rppp(window = Rectangle_window(), lambda = 100.0, nsim = 1, types = "type1", drop = TRUE)

  expect_equal(sample, other_sample)
})

test_that("Deduce number of types with only n.", {
  set.seed(42)
  sample <- rppp(lambda = c(1, 2, 3))
  set.seed(42)
  other_sample <- rppp(lambda = c(1, 2, 3), types = c("type1", "type2", "type3"))

  expect_equal(sample, other_sample)
})

test_that("Deduce number of types with only types.", {
  set.seed(42)
  sample <- rppp(types = c("a", "b", "c"))
  set.seed(42)
  other_sample <- rppp(lambda = c(100, 100, 100), types = c("a", "b", "c"))

  expect_equal(sample, other_sample)
})

test_that("Allow names to be given in n.", {
  set.seed(42)
  sample <- rppp(lambda = c(a = 1, b = 2, c = 3))
  set.seed(42)
  other_sample <- rppp(lambda = c(1, 2, 3), types = c("a", "b", "c"))

  expect_equal(sample, other_sample)
})

test_that("Allow n to be a list.", {
  set.seed(42)
  sample <- rppp(lambda = list(a = 1, b = 2, c = 3))
  set.seed(42)
  other_sample <- rppp(lambda = c(1, 2, 3), types = c("a", "b", "c"))

  expect_equal(sample, other_sample)
})

test_that("Allow types to be a list.", {
  set.seed(42)
  sample <- rppp(types = list("a", "b", "c"))
  set.seed(42)
  other_sample <- rppp(lambda = c(100, 100, 100), types = c("a", "b", "c"))

  expect_equal(sample, other_sample)
})

test_that("Look for inconsistencies in arguments n and types.", {
  expect_error(rppp(lambda = c(1, 1), types = c("a", "b", "c")))
  expect_error(rppp(lambda = c(1, 1), types = c("a")))

  expect_error(rppp(lambda = c(1, 1), types = list("a", "b", "c")))
  expect_error(rppp(lambda = c(1, 1), types = list("a")))

  expect_error(rppp(lambda = list(1, 1), types = c("a", "b", "c")))
  expect_error(rppp(lambda = list(1, 1), types = c("a")))
})


test_that("Points are in correct range.", {
  x_range <- c(-0.5, 4.5)
  y_range <- c(-0.8, -0.5)
  lambda <- 1.5
  number_tests <- 100
  window <- Rectangle_window(x_range, y_range)

  expect_true_condition = TRUE
  for(i in seq_along(number_tests)) {
    sample <- rppp(window = window, lambda = lambda, nsim = 1)
    x <- x_coordinates(sample)
    y <- y_coordinates(sample)
    expect_true_condition <- expect_true_condition && all(x >= x_range[1] & x <= x_range[2])
    expect_true_condition <- expect_true_condition && all(y >= y_range[1] & y <= y_range[2])
  }
  expect_true(expect_true_condition)
})

test_that("Number of points populated by correct number of calls to rpois.", {
  lambda <- 150
  arbitrary_seed <- 42
  default_range <- c(0, 1)
  window <- Rectangle_window(default_range, default_range)

  set.seed(arbitrary_seed)
  expected_values <- rpois(1, lambda)
  set.seed(arbitrary_seed)
  output <- rppp(window = window, lambda = lambda, nsim = 1)
  output_values <- length(x_coordinates(output))
  expect_equal(expected_values, output_values)
})
