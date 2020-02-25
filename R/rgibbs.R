
#' Sample a multivariate Gibbs point processes
#'
#' @param window Simulation window.
#' @param alpha Repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.
#' @param gamma Medium range repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.
#' @param lambda A vector representing the intensities of the point processes.
#' Default is a vector of same size as types, filled with ones.
#' @param covariates Covariates, with an empty list as a default.
#' @param beta Fitted coefficients related to covariates. Default is square matrix of zeroes of the same
#' number of rows/columns as the covariates.
#' @param short_range Symmetric matrix of short range interaction radii. Filled with 0.1 times the window radius by default.
#' @param medium_range Symmetric matrix of medium range interaction radii. Filled with 0.1 times the window radius by default.
#' @param long_range Symmetric matrix of long range interaction radii. Filled with 0.2 times the window radius by default.
#' @param saturation Saturation parameter of the point process. Default is 2.
#' @param steps Number of steps in the Metropolis algorithm. If `steps = 0`, uses the coupling from the past algorithm instead.
#' Default is 0.
#' @param nsim Number of samples to generate. Default is 1.
#' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
#' @param model String representing the model to use You can check the currently authorised models with a call to `show_short_range_models()`.
#' @param medium_range_model String representing the model to use You can check the currently authorised models with a call to `show_medium_range_models()`.
#' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration. Default is TRUE.
#' @param mark_range Range of additional marks to give to the points.
#' @export
rgibbs <- function(window = NULL,
                   alpha = NULL,
                   gamma = NULL,
                   lambda = NULL,
                   covariates = NULL,
                   beta = NULL,
                   short_range = NULL,
                   medium_range = NULL,
                   long_range = NULL,
                   saturation = 2,
                   steps = 0,
                   nsim = 1,
                   types = NULL,
                   model = "square_bump",
                   medium_range_model = "square_exponential",
                   drop = TRUE,
                   mark_range = c(1.0, 1.0)) {
  if(is.null(window)) {
    window <- Rectangle_window()
  }
  # Make covariates im objects with proper names.
  covariates <- coerce_to_named_im_objects(covariates, "unnamed_covariate", window)

  rgibbs_cpp(window, alpha, lambda, covariates, beta, gamma, short_range, medium_range, long_range, saturation, steps, nsim, types, model, medium_range_model, drop, mark_range)
}
