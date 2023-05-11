library(ppjsdm)

context("Configuration")

test_that("Configuration is an S3 class", {
  expect_is(ppjsdm::Configuration(0), "Configuration")
})

test_that("Configuration accessors", {
  x_coordinates <- stats::runif(n = 3)
  y_coordinates <- stats::runif(n = 3)
  types <- factor(c("a", "b", "c"))
  marks <- stats::runif(n = 3)
  configuration <- ppjsdm::Configuration(x = x_coordinates, y = y_coordinates, types = types, marks = marks)

  expect_true(is.vector(ppjsdm::x_coordinates(configuration)))
  expect_true(all.equal(ppjsdm::x_coordinates(configuration), x_coordinates))

  expect_true(is.vector(ppjsdm::y_coordinates(configuration)))
  expect_true(all.equal(ppjsdm::y_coordinates(configuration), y_coordinates))

  expect_true(is.factor(ppjsdm::types(configuration)))
  expect_true(all.equal(ppjsdm::types(configuration), types))

  expect_true(is.vector(ppjsdm::marks(configuration)))
  expect_true(all.equal(ppjsdm::marks(configuration), marks))
})

test_that("Configuration warning on duplicate points", {
  x_coordinates <- c(0, 1, 1, 2)
  y_coordinates <- c(2, 1, 1, 0)
  types <- factor(c("a", "b", "b", "c"))
  marks <- c(3, 2, 2, 1)

  expect_warning(configuration <- Configuration(x = x_coordinates, y = y_coordinates, types = types, marks = marks))
  expect_warning(Configuration(configuration))
})

test_that("Configuration constructor default arguments", {
  n <- 3
  x_default <- rep(1, n)
  y_default <- rep(1, n)
  types_default <- factor(rep("default", n))
  marks_default <- rep(1, n)

  # Expecting a warning about duplicate points here and in all tests below
  expect_warning(configuration_default <- ppjsdm::Configuration(x = x_default,
                                                                y = y_default,
                                                                types = types_default,
                                                                marks = marks_default))
  expect_warning(expect_equal(ppjsdm::Configuration(x = x_default), configuration_default))
  expect_warning(expect_equal(ppjsdm::Configuration(y = y_default), configuration_default))
  expect_warning(expect_equal(ppjsdm::Configuration(types = types_default), configuration_default))
  expect_warning(expect_equal(ppjsdm::Configuration(marks = marks_default), configuration_default))
})

test_that("Empty configuration", {
  expect_equal(ppjsdm::Configuration(x = vector(mode = "numeric", length = 0),
                                     y = vector(mode = "numeric", length = 0)),
               ppjsdm::Configuration())
})

test_that("Configuration copy-constructor", {
  x <- stats::runif(n = 1)
  y <- stats::runif(n = 1)
  types <- factor(c("a"))
  marks <- stats::runif(n = 1)

  configuration <- ppjsdm::Configuration(x = x, y = y, types = types, marks = marks)
  configuration_copy <- ppjsdm::Configuration(configuration)

  expect_identical(ppjsdm::x_coordinates(configuration_copy), x)
  expect_identical(ppjsdm::y_coordinates(configuration_copy), y)
  expect_identical(ppjsdm::types(configuration_copy), types)
  expect_identical(ppjsdm::marks(configuration_copy), marks)
})

test_that("Configuration constructor from matrix", {
  configuration <- ppjsdm::Configuration(cbind(c(0, 1), c(2, 3)))

  expect_identical(ppjsdm::x_coordinates(configuration), c(0, 1))
  expect_identical(ppjsdm::y_coordinates(configuration), c(2, 3))

  configuration <- ppjsdm::Configuration(cbind(c(0, 1), c(2, 3), c(4, 5)))

  expect_identical(ppjsdm::x_coordinates(configuration), c(0, 1))
  expect_identical(ppjsdm::y_coordinates(configuration), c(2, 3))
  expect_identical(ppjsdm::marks(configuration), c(4, 5))
})

test_that("Configuration constructor from data frame or list", {
  configuration <- ppjsdm::Configuration(data.frame(x = c(0, 1), y = c(2, 3)))

  expect_identical(ppjsdm::x_coordinates(configuration), c(0, 1))
  expect_identical(ppjsdm::y_coordinates(configuration), c(2, 3))

  configuration <- ppjsdm::Configuration(data.frame(x = c(0, 1), y = c(2, 3), marks = c(4, 5)))

  expect_identical(ppjsdm::x_coordinates(configuration), c(0, 1))
  expect_identical(ppjsdm::y_coordinates(configuration), c(2, 3))
  expect_identical(ppjsdm::marks(configuration), c(4, 5))

  configuration <- ppjsdm::Configuration(list(x = c(0, 1), y = c(2, 3), marks = c(4, 5), types = c(6, 7)))

  expect_identical(ppjsdm::x_coordinates(configuration), c(0, 1))
  expect_identical(ppjsdm::y_coordinates(configuration), c(2, 3))
  expect_identical(ppjsdm::marks(configuration), c(4, 5))
  expect_identical(ppjsdm::types(configuration), factor(c(6, 7)))
})

test_that("Configuration constructor from spatstat.geom::ppp", {
  x <- stats::runif(n = 2)
  y <- stats::runif(n = 2)
  types <- factor(c("a", "b"))
  configuration_spatstat <- ppjsdm::as.Configuration(spatstat.geom::ppp(x = x,
                                                                        y = y,
                                                                        marks = types,
                                                                        W = spatstat.geom::owin(c(0, 1), c(0, 3))))
  configuration <- ppjsdm::Configuration(x = x, y = y, types = types)
  expect_identical(configuration, configuration_spatstat)
})

test_that("Configuration unsupported constructors", {
  expect_error(ppjsdm::Configuration(c(0), c(1, 2)))
  expect_error(ppjsdm::Configuration(c(1, 2), 0))
  expect_error(ppjsdm::Configuration(c('a', 'b'), c(1, 2)))
  expect_error(ppjsdm::Configuration(x = c(0, 0.2), y = c(1, 2), types = factor("a")))
  expect_error(ppjsdm::Configuration(x = c(0, 0.2), y = c(1, 2), types = factor(c("a", "a", "a"))))
})

test_that("Coerce type to factor", {
  x_coordinates <- stats::runif(n = 3)
  y_coordinates <- stats::runif(n = 3)
  default_tag <- "a"
  expect_equal(ppjsdm::Configuration(x = x_coordinates, y = y_coordinates, types = default_tag),
               ppjsdm::Configuration(x = x_coordinates, y = y_coordinates, types = factor(c(default_tag, default_tag, default_tag))))
  default_tag <- 3.3
  expect_equal(ppjsdm::Configuration(x = x_coordinates, y = y_coordinates, types = default_tag),
               ppjsdm::Configuration(x = x_coordinates, y = y_coordinates, types = factor(c(default_tag, default_tag, default_tag))))
})

test_that("Configuration subscript operator", {
  x1 <- stats::runif(n = 2)
  x2 <- stats::runif(n = 2)
  x3 <- stats::runif(n = 2)
  x_coordinates <- c(x1[1], x2[1], x3[1])
  y_coordinates <- c(x1[2], x2[2], x3[2])
  types <- factor(c("a", "b", "a"))
  configuration <- ppjsdm::Configuration(x = x_coordinates, y = y_coordinates, types = types)

  expect_equal(configuration[1], ppjsdm::Configuration(x = c(x1[1], x3[1]), y = c(x1[2], x3[2]), types = factor(c("a", "a"))))
  expect_equal(configuration[2], ppjsdm::Configuration(x = c(x2[1]), y = c(x2[2]), types = factor(c("b"))))

  expect_equal(configuration[2], configuration[2.99])
  expect_equal(configuration[3], ppjsdm::Configuration(x = numeric(0), y = numeric(0), types = factor()))
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

test_that("Configuration number of points", {
  n1 <- 12
  n2 <- 41
  n3 <- 4
  n <- n1 + n2 + n3
  x_coordinates <- 0:(n - 1)
  y_coordinates <- 2:(n + 1)
  types <- factor(c(rep("a", n1), rep("b", n2), rep("c", n3)))
  configuration <- Configuration(x = x_coordinates, y = y_coordinates, types = types)

  expect_equal(length(configuration), n)
})

test_that("Convert to spatstat.geom::ppp", {
  x_coordinates <- stats::runif(n = 3)
  y_coordinates <- stats::runif(n = 3)
  types <- c("a", "b", "a")
  marks <- stats::runif(n = 3)
  expect_equal(as.ppp(ppjsdm::Configuration(x = x_coordinates, y = y_coordinates, types = types, marks = marks), W = Rectangle_window(c(0, 2), c(0, 4))),
               spatstat.geom::ppp(x = x_coordinates, y = y_coordinates, window = spatstat.geom::owin(c(0, 2), c(0, 4)), marks = factor(types)))
  types <- c(1, 2, 1)
  expect_equal(as.ppp(ppjsdm::Configuration(x = x_coordinates, y = y_coordinates, types = types, marks = marks), W = Rectangle_window(c(0, 2), c(0, 4))),
               spatstat.geom::ppp(x = x_coordinates, y = y_coordinates, window = spatstat.geom::owin(c(0, 2), c(0, 4)), marks = factor(types)))

})
