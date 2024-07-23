
#' Sample a birth-death point process which has as its stationary distribution a multivariate saturated pairwise interaction Gibbs point process.
#' This function has been tested well but its theoretical foundation is still work in progress. The parameters are also temporary and might have different names in later versions.
#'
#' @param window Observation window. Preferably a `ppjsdm` Window, such as `ppjsdm::Rectangle_window`, but also accepts `spatstat` `im` or `owin` objects.
#' @param nsim Number of samples to generate. Default is 1.
#' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
#' @param horizon Temporal horizon up to which we simulate
#' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration. Default is TRUE.
#' @param mark_range Range of additional marks to give to the points.
#' @param starting_configuration Optional configuration to start with when using the birth-death process.
#' @param nquad Number of quadrature points used to approximate the integral of the birth rate.
#' @param nthreads Number of threads to use. Default is 1.
#' @param birth_alpha Repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.
#' @param birth_gamma Medium range repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.
#' @param birth_beta0 A vector representing the log_intensities of the point processes.
#' Default is a vector of same size as types, filled with zeros
#' @param birth_covariates Covariates, with an empty list as a default.
#' @param birth_beta Fitted coefficients related to covariates. Default is square matrix of zeroes of the same
#' number of rows/columns as the covariates.
#' @param birth_short_range Symmetric matrix of short range interaction radii. Filled with 0.1 by default.
#' @param birth_medium_range Symmetric matrix of medium range interaction radii. Filled with 0 by default.
#' @param birth_long_range Symmetric matrix of long range interaction radii. Filled with 0 by default.
#' @param birth_saturation Saturation parameter of the point process. Default is 2.
#' @param birth_model String representing the model to use You can check the currently authorised models with a call to `show_short_range_models()`.
#' @param birth_medium_range_model String representing the model to use You can check the currently authorised models with a call to `show_medium_range_models()`.
#' @param death_alpha Repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.
#' @param death_gamma Medium range repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.
#' @param death_beta0 A vector representing the log_intensities of the point processes.
#' Default is a vector of same size as types, filled with zeros
#' @param death_covariates Covariates, with an empty list as a default.
#' @param death_beta Fitted coefficients related to covariates. Default is square matrix of zeroes of the same
#' number of rows/columns as the covariates.
#' @param death_short_range Symmetric matrix of short range interaction radii. Filled with 0.1 by default.
#' @param death_medium_range Symmetric matrix of medium range interaction radii. Filled with 0 by default.
#' @param death_long_range Symmetric matrix of long range interaction radii. Filled with 0 by default.
#' @param death_saturation Saturation parameter of the point process. Default is 2.
#' @param death_model String representing the model to use You can check the currently authorised models with a call to `show_short_range_models()`.
#' @param death_medium_range_model String representing the model to use You can check the currently authorised models with a call to `show_medium_range_models()`.
#' @importFrom spatstat.geom rQuasi
#' @export
#' @md
rbirth <- function(window,
                   nsim = 1,
                   types,
                   horizon = 1,
                   drop = TRUE,
                   mark_range = c(1.0, 1.0),
                   starting_configuration = NULL,
                   nquad = 1e3,
                   nthreads = 1,
                   birth_alpha,
                   birth_gamma,
                   birth_beta0,
                   birth_covariates,
                   birth_beta,
                   birth_short_range,
                   birth_medium_range,
                   birth_long_range,
                   birth_saturation,
                   birth_model,
                   birth_medium_range_model,
                   death_alpha,
                   death_gamma,
                   death_beta0,
                   death_covariates,
                   death_beta,
                   death_short_range,
                   death_medium_range,
                   death_long_range,
                   death_saturation,
                   death_model,
                   death_medium_range_model) {

  birth_parameters <- model_parameters(window = window,
                                       alpha = birth_alpha,
                                       gamma = birth_gamma,
                                       beta0 = birth_beta0,
                                       covariates = birth_covariates,
                                       beta = birth_beta,
                                       short_range = birth_short_range,
                                       medium_range = birth_medium_range,
                                       long_range = birth_long_range,
                                       saturation = birth_saturation,
                                       types = types,
                                       model = birth_model,
                                       medium_range_model = birth_medium_range_model)

  death_parameters <- model_parameters(window = window,
                                       alpha = death_alpha,
                                       gamma = death_gamma,
                                       beta0 = death_beta0,
                                       covariates = death_covariates,
                                       beta = death_beta,
                                       short_range = death_short_range,
                                       medium_range = death_medium_range,
                                       long_range = death_long_range,
                                       saturation = death_saturation,
                                       types = types,
                                       model = death_model,
                                       medium_range_model = death_medium_range_model)

  # TODO: Not sure about choosing dummy points with only one type. We end up computing Papangelou intensities
  # with these dummy points, should they not cover all possible types? Debugging other stuff right now,
  # so this will do for now.
  dummy <- as.Configuration(rQuasi(n = nquad, W = as.owin(birth_parameters$window), type = "Halton"))

  rbirth_cpp(nsim = nsim,
             horizon = horizon,
             types = birth_parameters$full_types,
             window = birth_parameters$window,
             drop = drop,
             nquad = nquad,
             seed = sample.int(.Machine$integer.max, 1),
             mark_range = mark_range,
             starting_configuration = starting_configuration,
             dummy = dummy,
             nthreads = nthreads,
             birth_alpha = birth_parameters$alpha,
             birth_beta0 = birth_parameters$beta0,
             birth_covariates = birth_parameters$covariates,
             birth_beta = birth_parameters$beta,
             birth_gamma = birth_parameters$gamma,
             birth_short_range = birth_parameters$short_range,
             birth_medium_range = birth_parameters$medium_range,
             birth_long_range = birth_parameters$long_range,
             birth_saturation = birth_parameters$saturation,
             birth_model = birth_parameters$model,
             birth_medium_range_model = birth_parameters$medium_range_model,
             death_alpha = death_parameters$alpha,
             death_beta0 = death_parameters$beta0,
             death_covariates = death_parameters$covariates,
             death_beta = death_parameters$beta,
             death_gamma = death_parameters$gamma,
             death_short_range = death_parameters$short_range,
             death_medium_range = death_parameters$medium_range,
             death_long_range = death_parameters$long_range,
             death_saturation = death_parameters$saturation,
             death_model = death_parameters$model,
             death_medium_range_model = death_parameters$medium_range_model)
}
