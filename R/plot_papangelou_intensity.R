#' Plot the Papangelou conditional intensity.
#'
#' @param window Observation window.
#' @param type Type of the point at which to evaluate the Papangelou conditional intensity.
#' @param configuration Configuration of points at which to evaluate the Papangelou conditional intensity.
#' @param model Model for short range interaction.
#' @param medium_range_model Model for medium range interaction.
#' @param alpha Short range repulsion parameter.
#' @param beta0 Log-intensity.
#' @param beta Covariates coefficients.
#' @param gamma Medium range repulsion parameter.
#' @param covariates Covariates.
#' @param short_range Matrix of short range distances.
#' @param medium_range Matrix of medium range distances.
#' @param long_range Matrix of long range distances.
#' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
#' @param saturation Saturation.
#' @param grid_steps Number of horizontal and vertical grid steps.
#' @param mark Mark of the point.
#' @param steps Nunber of steps in the Metropolis-Hastings simulation algorithm.
#' @param nthreads Maximum number of threads for parallel computing.
#' @param use_log Plot the logarithm of the Papangelou conditional intensity instead?
#' @importFrom spatstat.geom as.im as.owin as.ppp
#' @importFrom stats na.omit
#' @importFrom graphics plot
#' @export
plot_papangelou <- function(window,
                            type = 1,
                            mark = 1.0,
                            configuration,
                            model,
                            medium_range_model,
                            alpha,
                            beta0,
                            beta,
                            gamma,
                            covariates,
                            short_range,
                            medium_range,
                            long_range,
                            types,
                            saturation,
                            grid_steps = 1000,
                            steps = 0,
                            nthreads = 4,
                            use_log = FALSE) {

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

  if(missing(configuration)) {
    configuration <- rgibbs(window = parameters$window,
                            alpha = parameters$alpha,
                            beta0 = parameters$beta0,
                            beta = parameters$beta,
                            covariates = parameters$covariates,
                            short_range = parameters$short_range,
                            medium_range = parameters$medium_range,
                            long_range = parameters$long_range,
                            saturation = parameters$saturation,
                            model = parameters$model,
                            medium_range_model = parameters$medium_range_model,
                            gamma = parameters$gamma,
                            types = parameters$types,
                            steps = steps)
  }
  window <- as.owin(parameters$window)
  x_range <- window$xrange
  y_range <- window$yrange
  x_axis <- seq(from = x_range[1], to = x_range[2], length.out = grid_steps)
  y_axis <- seq(from = y_range[1], to = y_range[2], length.out = grid_steps)

  z <- outer(x_axis, y_axis, function(x, y) {
    compute_papangelou_cpp(x = x,
                           y = y,
                           type = rep(type, length(x)),
                           mark = rep(mark, length(x)),
                           configuration = as.Configuration(configuration),
                           model = parameters$model,
                           medium_range_model = parameters$medium_range_model,
                           alpha = parameters$alpha,
                           beta0 = parameters$beta0,
                           beta = parameters$beta,
                           gamma = parameters$gamma,
                           covariates = parameters$covariates,
                           short_range = parameters$short_range,
                           medium_range = parameters$medium_range,
                           long_range = parameters$long_range,
                           saturation = parameters$saturation,
                           nthreads = nthreads)
  })
  if(use_log) {
    z <- log(z)
    title <- "log-Papangelou conditional intensity"
  } else {
    title <- "Papangelou conditional intensity"
  }
  plot(as.im(t(z), W = window), main = title)
  plot(as.ppp(configuration, window), add = TRUE, cols = 'white')
}
