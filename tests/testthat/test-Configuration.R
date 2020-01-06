library(ppjsdm)

context("Configuration")

test_that("Configuration is an S3 class", {
  configuration <- Configuration(0, 1, factor(c("a")))

  expect_is(configuration, "Configuration")
})

test_that("Configuration accessors", {
  x_coordinates <- 0:2
  y_coordinates <- 3:5
  types <- factor(c("a", "b", "c"))
  configuration <- Configuration(x_coordinates, y_coordinates, types)

  expect_true(is.vector(x_coordinates(configuration)))
  expect_true(all.equal(x_coordinates(configuration), x_coordinates))

  expect_true(is.vector(y_coordinates(configuration)))
  expect_true(all.equal(y_coordinates(configuration), y_coordinates))

  expect_true(is.factor(types(configuration)))
  expect_true(all.equal(types(configuration), types))
})

test_that("Coerce type to factor", {
  x_coordinates <- 0:2
  y_coordinates <- 3:5
  default_tag <- "a"
  expect_equal(Configuration(x_coordinates, y_coordinates, default_tag),
               Configuration(x_coordinates, y_coordinates, factor(c(default_tag, default_tag, default_tag))))
  default_tag <- 3.3
  expect_equal(Configuration(x_coordinates, y_coordinates, default_tag),
               Configuration(x_coordinates, y_coordinates, factor(c(default_tag, default_tag, default_tag))))
})

test_that("Configuration copy-constructor", {
  configuration <- Configuration(0, 1, factor(c("a")))
  configuration_copy <- Configuration(configuration)

  expect_identical(x_coordinates(configuration_copy), 0)
  expect_identical(y_coordinates(configuration_copy), 1)
  expect_identical(types(configuration_copy), factor(c("a")))
})

test_that("Configuration constructor from matrix", {
  configuration <- Configuration(rbind(c(0, 1), c(2, 3)))

  expect_identical(x_coordinates(configuration), c(0, 2))
  expect_identical(y_coordinates(configuration), c(1, 3))
})

test_that("Configuration default type", {
  configuration <- Configuration(x = c(0, 1), y = c(2, 3))
  other_configuration <- Configuration(x = c(0, 1), y = c(2, 3), types = factor(c("default", "default")))

  expect_equal(configuration, other_configuration)
})

test_that("Configuration warning on duplicate points", {
  x_coordinates <- c(0, 1, 1, 2)
  y_coordinates <- c(2, 1, 1, 0)
  types <- factor(c("a", "b", "b", "c"))

  expect_warning(configuration <- Configuration(x = x_coordinates, y = y_coordinates, types = types))
  expect_warning(Configuration(configuration))
})

test_that("Configuration unsupported constructors", {
  expect_error(Configuration())
  expect_error(Configuration(c(0), c(1, 2)))
  expect_error(Configuration(c(1, 2), 0))
  expect_error(Configuration(c('a', 'b'), c(1, 2)))
  expect_error(Configuration(x = c(0, 0.2), y = c(1, 2), types = factor("a")))
  expect_error(Configuration(x = c(0, 0.2), y = c(1, 2), types = factor(c("a", "a", "a"))))
})


test_that("Configuration constructor specified arguments", {
  x_coordinates <- c(0, 1)
  y_coordinates <- c(2, 3)
  types <- factor(c("a", "b"))
  configuration <- Configuration(x = x_coordinates, y = y_coordinates, types = types)

  expect_identical(x_coordinates(configuration), x_coordinates)
  expect_identical(y_coordinates(configuration), y_coordinates)
  expect_identical(types(configuration), types)
})

test_that("Configuration constructor unspecified arguments", {
  x_coordinates <- c(0, 1)
  y_coordinates <- c(2, 3)
  types <- factor(c("a", "b"))
  configuration <- Configuration(x_coordinates, y_coordinates, types)

  expect_identical(x_coordinates(configuration), x_coordinates)
  expect_identical(y_coordinates(configuration), y_coordinates)
  expect_identical(types(configuration), types)
})

test_that("Configuration subscript operator", {
  x_coordinates <- 0:2
  y_coordinates <- 3:5
  types <- factor(c("a", "b", "a"))
  configuration <- Configuration(x = x_coordinates, y = y_coordinates, types = types)

  expect_equal(configuration[1], Configuration(x = c(0, 2), y = c(3, 5), types = factor(c("a", "a"))))
  expect_equal(configuration[2], Configuration(x = c(1), y = c(4), types = factor(c("b"))))

  expect_equal(configuration[2], configuration[2.99])
  expect_equal(configuration[3], Configuration(x = numeric(0), y = numeric(0), types = factor()))
})

test_that("Configuration number of points", {
  n1 <- 12
  n2 <- 41
  n3 <- 4
  n <- n1 + n2 + n3
  x_coordinates <- 0:(n - 1)
  y_coordinates <- 2:(n + 1)
  types <- factor(c(rep("a", n1), rep("b", n2), rep("c", n3)))
  configuration <- Configuration(x = x_coordinates, y = y_coordinates, types = types)

  expect_equal(get_number_points(configuration), list(a = n1, b = n2, c = n3))
  expect_equal(get_number_points(configuration, total = FALSE), list(a = n1, b = n2, c = n3))
  expect_equal(get_number_points(configuration, total = TRUE), n)
})

test_that("Configuration vectorised subscript operator", {
  x_coordinates <- 0:3
  y_coordinates <- 4:7
  types <- factor(c("a", "b", "a", "c"))
  configuration <- Configuration(x = x_coordinates, y = y_coordinates, types = types)

  expect_equal(configuration[1:2], Configuration(x = c(0, 1, 2), y = c(4, 5, 6), types = factor(c("a", "b", "a"))))
  expect_equal(configuration[2:3], Configuration(x = c(1, 3), y = c(5, 7), types = factor(c("b", "c"))))
  expect_equal(configuration[1:3], configuration)
  expect_equal(configuration[-1], configuration[2:3])
  expect_equal(configuration[c(1, 3)], Configuration(x = c(0, 2, 3), y = c(4, 6, 7), types = factor(c("a", "a", "c"))))
})
