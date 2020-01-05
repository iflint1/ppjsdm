library(ppjsdm)

context("Disk_window")

test_that("Disk_window is an S4 class", {
  window <- Disk_window()

  expect_is(window, "Disk_window")
})

test_that("Disk_window accessors", {
  centre <- c(-1, 2)
  radius <- 0.5
  window <- Disk_window(centre, radius)

  expect_true(is.vector(centre(window)))
  expect_true(all.equal(centre(window), centre))

  expect_true(is.numeric(radius(window)))
  expect_true(all.equal(radius(window), radius))
})


test_that("Disk_window default constructor", {
  default_centre <- c(0, 0)
  default_radius <- 1
  default_window <- Disk_window()

  expect_identical(centre(default_window), default_centre)
  expect_identical(radius(default_window), default_radius)
})

test_that("Disk_window specified arguments constructor", {
  centre <- c(-1, 2)
  radius <- 0.5
  window_specified_arguments <- Disk_window(radius = radius, centre = centre)

  expect_identical(centre(window_specified_arguments), centre)
  expect_identical(radius(window_specified_arguments), radius)
})

test_that("Disk_window unspecified arguments constructors", {
  centre <- c(-1, 2)
  radius <- 1.5
  window_unspecified_arguments <- Disk_window(centre, radius)

  expect_identical(centre(window_unspecified_arguments), centre)
  expect_identical(radius(window_unspecified_arguments), radius)
})

test_that("Disk_window one argument constructor", {
  default_radius <- 1
  centre <- c(-1, 2)
  window_one_argument <- Disk_window(centre)

  expect_identical(centre(window_one_argument), centre)
  expect_identical(radius(window_one_argument), default_radius)
})

test_that("Disk_window copy constructor", {
  centre <- c(-1, 2)
  radius <- 0.5
  window_specified_arguments <- Disk_window(centre = centre, radius = radius)
  copied_window <- Disk_window(window_specified_arguments)

  expect_identical(centre(window_specified_arguments), centre(copied_window))
  expect_identical(radius(window_specified_arguments), radius(copied_window))
})

test_that("Disk_window unsupported constructors", {
  default_centre <- c(0, 0)

  expect_error(Disk_window(c(0, 1, 3)))
  expect_error(Disk_window(c('a', 'b')))
  expect_error(Disk_window(7))

  expect_error(Disk_window(default_centre, c(0, 1, 3)))
  expect_error(Disk_window(default_centre, c('a', 'b')))
  expect_error(Disk_window(default_centre, c(0, 0)))
  expect_error(Disk_window(default_centre, -1))
})
