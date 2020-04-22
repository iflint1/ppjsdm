#' Plot the Papangelou conditional intensity.
#'
#' @param window Observation window.
#' @param type Type of the point at which to evaluate the Papangelou conditional intensity.
#' @param configuration Configuration of points at which to evaluate the Papangelou conditional intensity.
#' @param model Model for short range interaction.
#' @param medium_range_model Model for medium range interaction.
#' @param alpha Short range repulsion parameter.
#' @param lambda Intensity.
#' @param beta Covariates coefficients.
#' @param gamma Medium range repulsion parameter.
#' @param covariates Covariates.
#' @param short_range Matrix of short range distances.
#' @param medium_range Matrix of medium range distances.
#' @param long_range Matrix of long range distances.
#' @param saturation Saturation.
#' @param grid_steps Number of horizontal and vertical grid steps.
#' @param mark Mark of the point.
#' @param steps Nunber of steps in the Metropolis-Hastings simulation algorithm.
#' @importFrom spatstat as.im as.owin as.ppp
#' @export
plot_papangelou <- function(window, type, configuration, model, medium_range_model, alpha, lambda,
                            beta, gamma, covariates, short_range, medium_range, long_range, saturation, grid_steps = 1000, mark = 1.0, steps = 0) {
  if(missing(configuration)) {
    configuration <- rgibbs(window = window,
                            alpha = alpha,
                            lambda = lambda,
                            beta = beta,
                            covariates = covariates,
                            short_range = short_range,
                            medium_range = medium_range,
                            long_range = long_range,
                            saturation = saturation,
                            model = model,
                            medium_range_model = medium_range_model,
                            gamma = gamma,
                            steps = steps)
  }
  window <- as.owin(window)
  x_range <- window$xrange
  y_range <- window$yrange
  x_axis <- seq(from = x_range[1], to = x_range[2], length.out = grid_steps)
  y_axis <- seq(from = y_range[1], to = y_range[2], length.out = grid_steps)
  z <- outer(x_axis, y_axis, function(x, y) {
    compute_papangelou(configuration, x, y, type, model, medium_range_model, alpha,
                       lambda, beta, gamma, covariates, short_range, medium_range, long_range, saturation, mark)
  })
  plot(as.im(t(z), W = window))
  plot(as.ppp(configuration, window), add = TRUE)
}
