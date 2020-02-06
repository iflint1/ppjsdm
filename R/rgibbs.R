
#' Sample a multivariate Gibbs point processes
#'
#' @param window Simulation window.
#' @param alpha Repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.
#' @param lambda A vector representing the intensities of the point processes.
#' Default is a vector of same size as types, filled with ones.
#' @param covariates Covariates, with an empty list as a default.
#' @param beta Fitted coefficients related to covariates. Default is square matrix of zeroes of the same
#' number of rows/columns as the covariates.
#' @param radius Symmetric matrix of interaction radii. Filled with 0.1 times the window radius by default.
#' @param saturation Saturation parameter of the point process. Default is 2.
#' @param steps Number of steps in the Metropolis algorithm. If `steps = 0`, uses the coupling from the past algorithm instead.
#' Default is 0.
#' @param nsim Number of samples to generate. Default is 1.
#' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
#' @param model String representing the model to simulate from. You can check the currently authorised models with a call to `show_model()`.
#' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration. Default is TRUE.
#' @export
rgibbs <- function(window = NULL,
                   alpha = NULL,
                   lambda = NULL,
                   covariates = NULL,
                   beta = NULL,
                   radius = NULL,
                   saturation = 2,
                   steps = 0,
                   nsim = 1,
                   types = NULL,
                   model = "square_bump",
                   drop = TRUE) {
  if(is.null(window)) {
    window <- Rectangle_window()
  }
  # Make covariates im objects with proper names.
  covariates <- coerce_to_named_im_objects(covariates, "unnamed_covariate", window)

  rgibbs_cpp(window, alpha, lambda, covariates, beta, radius, saturation, steps, nsim, types, model, drop)
}
