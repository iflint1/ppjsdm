#' Sample a multivariate saturated pairwise interaction Gibbs point processes.
#'
#' IMPORTANT: Check ?rgibbs.default for the documentation.
#'
#' @param ... Parameters to be forwarded to the relevant method. Currently, either
#' the parameters of the Gibbs point process, or a fit object obtained from running `gibbsm`.
#' @export
rgibbs <- function(...) {
  UseMethod("rgibbs")
}

#' Sample a multivariate saturated pairwise interaction Gibbs point processes with given parameters.
#'
#' @param window Observation window. Preferably a `ppjsdm` Window, such as `ppjsdm::Rectangle_window`, but also accepts `spatstat` `im` or `owin` objects.
#' @param alpha Repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.
#' @param gamma Medium range repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.
#' @param beta0 A vector representing the log_intensities of the point processes.
#' Default is a vector of same size as types, filled with zeros
#' @param covariates Covariates, with an empty list as a default.
#' @param beta Fitted coefficients related to covariates. Default is square matrix of zeroes of the same
#' number of rows/columns as the covariates.
#' @param short_range Symmetric matrix of short range interaction radii. Filled with 0.1 by default.
#' @param medium_range Symmetric matrix of medium range interaction radii. Filled with 0 by default.
#' @param long_range Symmetric matrix of long range interaction radii. Filled with 0 by default.
#' @param saturation Saturation parameter of the point process. Default is 2.
#' @param steps Number of steps in the Metropolis algorithm. If `steps = 0`, uses the coupling from the past algorithm instead.
#' Default is 0.
#' @param nsim Number of samples to generate. Default is 1.
#' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
#' @param model String representing the model to use You can check the currently authorised models with a call to `show_short_range_models()`.
#' @param medium_range_model String representing the model to use You can check the currently authorised models with a call to `show_medium_range_models()`.
#' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration. Default is TRUE.
#' @param mark_range Range of additional marks to give to the points.
#' @param only_simulate_these_types Parameter used to do conditional simulation, only drawing from some of the types (usually conditional upon the others).
#' @param conditional_configuration Simulate conditional on this configuration.
#' @param starting_configuration Optional configuration to start with when using the Metropolis-Hastings algorithm (steps > 0).
#' @param nthreads Number of threads to use. Default is 1.
#' @param debug Output debugging information (current only used with Metropolis-Hastings)?
#' @param ... Not used.
#' @export
#' @method rgibbs default
#' @md
rgibbs.default <- function(window,
                           alpha,
                           gamma,
                           beta0,
                           covariates,
                           beta,
                           short_range,
                           medium_range,
                           long_range,
                           saturation,
                           types,
                           model,
                           medium_range_model,
                           steps = 0,
                           nsim = 1,
                           drop = TRUE,
                           mark_range = c(1.0, 1.0),
                           only_simulate_these_types = NULL,
                           conditional_configuration = NULL,
                           starting_configuration = NULL,
                           nthreads = 1,
                           debug = FALSE,
                           ...) {
  parameters <- model_parameters(window = window,
                                 alpha = alpha,
                                 gamma = gamma,
                                 beta0 = beta0,
                                 covariates = covariates,
                                 beta = beta,
                                 short_range = short_range,
                                 medium_range = medium_range,
                                 long_range = long_range,
                                 saturation = saturation,
                                 types = types,
                                 model = model,
                                 medium_range_model = medium_range_model)

  if(is.null(only_simulate_these_types)) {
    only_simulate_these_types <- seq_len(length(parameters$beta0))
  }
  only_simulate_these_types <- only_simulate_these_types - 1

  if(!is.null(starting_configuration)) {
    starting_configuration <- as.Configuration(starting_configuration)
  }

  rgibbs_cpp(window = parameters$window,
             alpha = parameters$alpha,
             beta0 = parameters$beta0,
             covariates = parameters$covariates,
             beta = parameters$beta,
             gamma = parameters$gamma,
             short_range = parameters$short_range,
             medium_range = parameters$medium_range,
             long_range = parameters$long_range,
             saturation = parameters$saturation,
             steps = steps,
             nsim = nsim,
             types = parameters$types,
             model = parameters$model,
             medium_range_model = parameters$medium_range_model,
             drop = drop,
             mark_range = mark_range,
             only_simulate_these_types = only_simulate_these_types,
             conditional_configuration = conditional_configuration,
             starting_configuration = starting_configuration,
             seed = sample.int(.Machine$integer.max, 1),
             nthreads = nthreads,
             debug = debug)
}

#' Sample a multivariate Gibbs point processes from a fit object.
#'
#' @param fit Fit object obtained by running gibbsm.
#' @param window Observation window. Preferably a `ppjsdm` Window, such as `ppjsdm::Rectangle_window`, but also accepts `spatstat` `im` or `owin` objects.
#' @param alpha Repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.
#' @param gamma Medium range repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.
#' @param beta0 A vector representing the log_intensities of the point processes.
#' Default is a vector of same size as types, filled with zeros
#' @param covariates Covariates, with an empty list as a default.
#' @param beta Fitted coefficients related to covariates. Default is square matrix of zeroes of the same
#' number of rows/columns as the covariates.
#' @param short_range Symmetric matrix of short range interaction radii. Filled with 0.1 by default.
#' @param medium_range Symmetric matrix of medium range interaction radii. Filled with 0 by default.
#' @param long_range Symmetric matrix of long range interaction radii. Filled with 0 by default.
#' @param saturation Saturation parameter of the point process. Default is 2.
#' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
#' @param model String representing the model to use You can check the currently authorised models with a call to `show_short_range_models()`.
#' @param medium_range_model String representing the model to use You can check the currently authorised models with a call to `show_medium_range_models()`.
#' @param mark_range Range of additional marks to give to the points.
#' @param ... Forwarded to rgibbs.
#' @export
#' @method rgibbs gibbsm
#' @md
rgibbs.gibbsm <- function(fit,
                          window = fit$window,
                          alpha = fit$coefficients$alpha,
                          gamma = fit$coefficients$gamma,
                          beta0 = fit$coefficients$beta0,
                          covariates = fit$parameters$covariates,
                          beta = fit$coefficients$beta,
                          short_range = fit$coefficients$short_range,
                          medium_range = fit$coefficients$medium_range,
                          long_range = fit$coefficients$long_range,
                          saturation = fit$parameters$saturation,
                          types = fit$type_names,
                          model = fit$parameters$model,
                          medium_range_model = fit$parameters$medium_range_model,
                          mark_range = fit$mark_range,
                          ...) {
  rgibbs(window = window,
         alpha = alpha,
         gamma = gamma,
         beta0 = beta0,
         covariates = covariates,
         beta = beta,
         short_range = short_range,
         medium_range = medium_range,
         long_range = long_range,
         saturation = saturation,
         types = types,
         model = model,
         medium_range_model = medium_range_model,
         mark_range = mark_range,
         ...)
}
