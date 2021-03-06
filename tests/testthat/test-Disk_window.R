library(ppjsdm)

context("Disk_window")

test_that("Disk_window has correct type", {
  expect_is(Disk_window(), "Disk_window")
})


test_that("Disk_window default constructor", {
  default_centre <- c(0, 0)
  default_radius <- 1
  default_window <- ppjsdm::Disk_window()
  window <- ppjsdm::Disk_window(centre = default_centre, radius = default_radius)

  expect_identical(default_window, window)
})

test_that("Disk_window one argument constructor", {
  default_radius <- 1
  centre <- stats::runif(n = 2)
  window_one_argument <- ppjsdm::Disk_window(centre)
  window <- ppjsdm::Disk_window(centre = centre, radius = default_radius)

  expect_identical(window_one_argument, window)
})

test_that("Disk_window copy constructor", {
  centre <- stats::runif(n = 2)
  radius <- 0.5
  window_specified_arguments <- ppjsdm::Disk_window(centre = centre, radius = radius)
  copied_window <- ppjsdm::Disk_window(window_specified_arguments)

  expect_identical(window_specified_arguments, copied_window)
})

test_that("Disk_window unsupported constructors", {
  default_centre <- c(0, 0)

  expect_error(ppjsdm::Disk_window(c(0, 1, 3)))
  expect_error(ppjsdm::Disk_window(c('a', 'b')))
  expect_error(ppjsdm::Disk_window(7))

  expect_error(ppjsdm::Disk_window(default_centre, c(0, 1, 3)))
  expect_error(ppjsdm::Disk_window(default_centre, c('a', 'b')))
  expect_error(ppjsdm::Disk_window(default_centre, c(0, 0)))
  expect_error(ppjsdm::Disk_window(default_centre, -1))
})
