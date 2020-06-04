library(ppjsdm)
library(methods)

context("rbinomialpp")

test_that("Return type.", {
  n <- 19
  nsim <- 14
  sample <- rbinomialpp(n = n, nsim = nsim, drop = FALSE)

  expect_true(is.list(sample))
  expect_true(methods::is(obj = sample[[1]], class2 = "Configuration"))
  expect_equal(length(sample), nsim)
  expect_equal(length(x_coordinates(sample[[1]])), n)

  sample <- rbinomialpp(n = n, nsim = 1, drop = TRUE)

  expect_true(methods::is(obj = sample, class2 = "Configuration"))
  expect_equal(length(x_coordinates(sample)), n)
})

test_that("rppp wrong window type.", {
  window <- structure(list(x_range = x_range, y_range = y_range))
  expect_error(rbinomialpp(window = window))
})

test_that("Default arguments.", {
  set.seed(42)
  sample <- rbinomialpp()
  set.seed(42)
  other_sample <- rbinomialpp(window = Rectangle_window(), n = 1, nsim = 1, types = "type1", drop = TRUE)

  expect_equal(sample, other_sample)
})

test_that("Deduce number of types with only n.", {
  set.seed(42)
  sample <- rbinomialpp(n = c(1, 2, 3))
  set.seed(42)
  other_sample <- rbinomialpp(n = c(1, 2, 3), types = c("type1", "type2", "type3"))

  expect_equal(sample, other_sample)
})

test_that("Deduce number of types with only types.", {
  set.seed(42)
  sample <- rbinomialpp(types = c("a", "b", "c"))
  set.seed(42)
  other_sample <- rbinomialpp(n = c(1, 1, 1), types = c("a", "b", "c"))

  expect_equal(sample, other_sample)
})

test_that("Allow names to be given in n.", {
  set.seed(42)
  sample <- rbinomialpp(n = c(a = 1, b = 2, c = 3))
  set.seed(42)
  other_sample <- rbinomialpp(n = c(1, 2, 3), types = c("a", "b", "c"))

  expect_equal(sample, other_sample)
})

test_that("Allow n to be a list.", {
  set.seed(42)
  sample <- rbinomialpp(n = list(a = 1, b = 2, c = 3))
  set.seed(42)
  other_sample <- rbinomialpp(n = c(1, 2, 3), types = c("a", "b", "c"))

  expect_equal(sample, other_sample)
})

test_that("Allow types to be a list.", {
  set.seed(42)
  sample <- rbinomialpp(types = list("a", "b", "c"))
  set.seed(42)
  other_sample <- rbinomialpp(n = c(1, 1, 1), types = c("a", "b", "c"))

  expect_equal(sample, other_sample)
})

test_that("Look for inconsistencies in arguments n and types.", {
  expect_error(rbinomialpp(n = c(1, 1), types = c("a", "b", "c")))
  expect_error(rbinomialpp(n = c(1, 1), types = c("a")))

  expect_error(rbinomialpp(n = c(1, 1), types = list("a", "b", "c")))
  expect_error(rbinomialpp(n = c(1, 1), types = list("a")))

  expect_error(rbinomialpp(n = list(1, 1), types = c("a", "b", "c")))
  expect_error(rbinomialpp(n = list(1, 1), types = c("a")))
})

test_that("Is in correct range (rectangle window).", {
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

test_that("Is in correct range (disk window).", {
  centre <- c(-0.5, 4.5)
  radius <- 1.9
  number_tests <- 100
  number_samples <- 10
  window <- Disk_window(centre = centre, radius = radius)

  expect_true_condition = TRUE
  for(i in seq_along(number_tests)) {
    sample <- rbinomialpp(window = window, n = number_samples, nsim = 1)
    x <- x_coordinates(sample)
    y <- y_coordinates(sample)
    expect_true_condition <- expect_true_condition && all((x - centre[1]) * (x - centre[1]) + (y - centre[2]) * (y - centre[2]) <= radius * radius)
  }
  expect_true(expect_true_condition)
})

test_that("Has correct number of points.", {
  x_range <- c(-0.5, 4.5)
  y_range <- c(-0.8, -0.5)
  window <- Rectangle_window(x_range, y_range)
  n <- c(7, 2, 9, 4)

  sample <- rbinomialpp(window = window, n = n, nsim = 1)

  ltypes <- levels(types(sample))
  num <- sapply(ltypes, function(l) length(x_coordinates(sample)[types(sample) == l]), USE.NAMES = FALSE)

  expect_equal(num, c(n[1], n[2], n[3], n[4]))
  expect_equal(length(sample[1]), n[1])
  expect_equal(length(sample[2]), n[2])
  expect_equal(length(sample[3]), n[3])
  expect_equal(length(sample[4]), n[4])
})

test_that("Samples from the window populated by right amount of calls to runif.", {
  default_range <- c(0, 1)
  arbitrary_seed <- 42
  window <- Rectangle_window(default_range, default_range)

  set.seed(arbitrary_seed)
  expected_values <- runif(n = 6, min = 0, max = 1) # Note: need two additional uniforms to sample the marks.
  set.seed(arbitrary_seed)
  output <- rbinomialpp(window = window, n = 2, nsim = 1)
  output_values <- c(x_coordinates(output), y_coordinates(output))
  expect_true(all(output_values %in% expected_values))
})
