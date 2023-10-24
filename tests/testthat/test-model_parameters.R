library(ppjsdm)

context("Model default parameters")

test_that("Default model parameters.", {
  no_arguments <- ppjsdm::model_parameters()

  all_arguments <- ppjsdm::model_parameters(window = ppjsdm::Rectangle_window(),
                                            alpha = 0.,
                                            gamma = 0.,
                                            beta0 = 0.,
                                            covariates = list(),
                                            beta = matrix(0, ncol = 0, nrow = 1),
                                            short_range = 0.1,
                                            medium_range = 0.,
                                            long_range = 0.,
                                            saturation = 2,
                                            types = "type1",
                                            model = "square_bump",
                                            medium_range_model = "square_exponential",
                                            default_number_types = 1)

  expect_identical(no_arguments, all_arguments)
})

test_that("Correctly detect the number of potentials with 1 type", {
  explicit <- ppjsdm::model_parameters(window = ppjsdm::Rectangle_window(),
                                       alpha = list(0., 0.),
                                       gamma = 0.,
                                       beta0 = 0.,
                                       covariates = list(),
                                       beta = matrix(0, ncol = 0, nrow = 1),
                                       short_range = list(0.05, 0.1),
                                       medium_range = 0.,
                                       long_range = 0.,
                                       saturation = 2,
                                       types = "type1",
                                       model = list("square_bump", "square_bump"),
                                       medium_range_model = "square_exponential",
                                       default_number_types = 1)

  implicit <- ppjsdm::model_parameters(alpha = list(0., 0.))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(alpha = c(0., 0.))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(short_range = list(0.05, 0.1))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(short_range = c(0.05, 0.1))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(model = list("square_bump", "square_bump"))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(model = c("square_bump", "square_bump"))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(alpha = c(0., 0.), model = "square_bump")
  expect_identical(implicit, explicit)
})

test_that("Correctly detect the number of potentials with 3 types", {
  explicit <- ppjsdm::model_parameters(window = ppjsdm::Rectangle_window(),
                                       alpha = list(matrix(0., 3, 3), matrix(0., 3, 3)),
                                       gamma = matrix(0., 3, 3),
                                       beta0 = c(0., 0., 0.),
                                       covariates = list(),
                                       beta = matrix(0, ncol = 0, nrow = 3),
                                       short_range = list(matrix(0.05, 3, 3), matrix(0.1, 3, 3)),
                                       medium_range = matrix(0., 3, 3),
                                       long_range = matrix(0., 3, 3),
                                       saturation = 2,
                                       types = c("type1", "type2", "type3"),
                                       model = list("square_bump", "square_bump"),
                                       medium_range_model = "square_exponential",
                                       default_number_types = 1)

  implicit <- ppjsdm::model_parameters(alpha = list(matrix(0., 3, 3), matrix(0., 3, 3)))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(alpha = list(matrix(0., 3, 3), 0.))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(alpha = list(0., matrix(0., 3, 3)))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(short_range = list(matrix(0.05, 3, 3), matrix(0.1, 3, 3)))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(short_range = list(matrix(0.05, 3, 3), 0.1))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(short_range = list(0.05, matrix(0.1, 3, 3)))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(short_range = list(matrix(0.05, 3, 3), matrix(0.1, 3, 3)),
                                       model = "square_bump")
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(short_range = list(matrix(0.05, 3, 3), matrix(0.1, 3, 3)),
                                       alpha = 0.)
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(alpha = c(0., 0.),
                                       gamma = matrix(0., 3, 3))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(short_range = c(0.05, 0.1),
                                       gamma = matrix(0., 3, 3))
  expect_identical(implicit, explicit)
})

test_that("Figure out the number of types indirectly", {
  explicit <- ppjsdm::model_parameters(window = ppjsdm::Rectangle_window(),
                                       alpha = matrix(0., 3, 3),
                                       gamma = matrix(0., 3, 3),
                                       beta0 = c(0., 0., 0.),
                                       covariates = list(),
                                       beta = matrix(0, ncol = 0, nrow = 3),
                                       short_range = matrix(0.1, 3, 3),
                                       medium_range = matrix(0., 3, 3),
                                       long_range = matrix(0., 3, 3),
                                       saturation = 2,
                                       types = c("type1", "type2", "type3"),
                                       model = "square_bump",
                                       medium_range_model = "square_exponential",
                                       default_number_types = 1)

  implicit <- ppjsdm::model_parameters(types = c("type1", "type2", "type3"))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(alpha = matrix(0., 3, 3))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(gamma = matrix(0., 3, 3))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(beta0 = c(0., 0., 0.))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(short_range = matrix(0.1, 3, 3))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(medium_range = matrix(0., 3, 3))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(long_range = matrix(0., 3, 3))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(something = c(1, 2, 3))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(something = matrix(0, 3, 3))
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(default_number_types = 3)
  expect_identical(implicit, explicit)

  implicit <- ppjsdm::model_parameters(covariates = list(temp = function(x, y) x),
                                       beta = matrix(0, ncol = 1, nrow = 3))
  explicit <- ppjsdm::model_parameters(window = ppjsdm::Rectangle_window(),
                                       alpha = matrix(0., 3, 3),
                                       gamma = matrix(0., 3, 3),
                                       beta0 = c(0., 0., 0.),
                                       covariates = list(temp = function(x, y) x),
                                       beta = matrix(0, ncol = 1, nrow = 3),
                                       short_range = matrix(0.1, 3, 3),
                                       medium_range = matrix(0., 3, 3),
                                       long_range = matrix(0., 3, 3),
                                       saturation = 2,
                                       types = c("type1", "type2", "type3"),
                                       model = "square_bump",
                                       medium_range_model = "square_exponential",
                                       default_number_types = 1)
  expect_identical(implicit, explicit)
})
