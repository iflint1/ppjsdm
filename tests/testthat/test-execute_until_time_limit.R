library(ppjsdm)

context("Execute until time limit")

test_that("Returns initial objects when time limit is large", {
  func <- function(x) {
    Sys.sleep(0.1)
    x
  }

  objects <- lapply(1:4, function(i) i)
  expect_equal(objects, ppjsdm:::execute_until_time_limit(objects = objects, func = func, time_limit = 1, unit = "secs"))

  objects <- lapply(c(1.2, 1.5, 1.5), function(i) i)
  expect_equal(objects, ppjsdm:::execute_until_time_limit(objects = objects, func = func, time_limit = Inf, unit = "secs"))

  objects <- lapply(c(TRUE, FALSE), function(i) i)
  expect_equal(objects, ppjsdm:::execute_until_time_limit(objects = objects, func = func, time_limit = 0.5, unit = "secs"))

  objects <- lapply(c("a", "b", "c"), function(i) i)
  expect_equal(objects, ppjsdm:::execute_until_time_limit(objects = objects, func = func, time_limit = 1, unit = "mins"))
})

test_that("Returns first object when time limit is small", {
  func <- function(x) {
    Sys.sleep(0.1)
    x
  }

  objects <- lapply(1:4, function(i) i)
  expect_equal(objects[1], ppjsdm:::execute_until_time_limit(objects = objects, func = func, time_limit = 0, unit = "secs"))

  objects <- lapply(c(1.2, 1.5, 1.5), function(i) i)
  expect_equal(objects[1], ppjsdm:::execute_until_time_limit(objects = objects, func = func, time_limit = 0.01, unit = "secs"))
})

test_that("Stays within the specified time limit", {
  delta <- 0.001
  func <- function(x) {
    Sys.sleep(delta)
    x
  }

  objects <- lapply(1:100, function(i) i)
  for(i in seq_len(100)) {
    ret <- ppjsdm:::execute_until_time_limit(objects = objects, func = func, time_limit = i * delta, unit = "secs")
    expect(length(ret) <= i, failure_message = paste0("Too many function executions, i = ", i, "."))
  }
})
