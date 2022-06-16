#' Plot the Papangelou conditional intensity.
#'
#' @param ... Parameters to be forwarded to the relevant method. Currently, either
#' the parameters of the Gibbs point process, or a fit object obtained from running `gibbsm`.
#' @export
plot_papangelou <- function(...) {
  UseMethod("plot_papangelou")
}

#' Plot the Papangelou conditional intensity with given parameters.
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
#' @param grid_steps Vector of length 2, representing the number of horizontal and vertical grid steps.
#' @param mark Mark of the point.
#' @param steps Nunber of steps in the Metropolis-Hastings simulation algorithm.
#' @param nthreads Maximum number of threads for parallel computing.
#' @param use_log Plot the logarithm of the Papangelou conditional intensity instead?
#' @param use_ggplot Use ggplot for fancier plot?
#' @param return_papangelou Should we return the Papangelou intensity itself instead of plotting it?
#' @param limits Limits for values of the Papamgelou conditional intensity in plotting functions.
#' @param ... Ignored.
#' @importFrom colorspace scale_fill_continuous_sequential
#' @importFrom ggplot2 aes coord_equal element_text geom_tile ggplot ggtitle guide_legend guides labs scale_color_manual scale_shape_manual scale_x_continuous scale_y_continuous theme theme_minimal xlab ylab
#' @importFrom graphics plot
#' @importFrom spatstat.geom as.im as.owin as.ppp
#' @importFrom stats na.omit
#' @export
#' @method plot_papangelou default
plot_papangelou.default <- function(window,
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
                                    grid_steps = c(1000, 1000),
                                    steps = 0,
                                    nthreads = 4,
                                    use_log = FALSE,
                                    use_ggplot = TRUE,
                                    return_papangelou = FALSE,
                                    limits,
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
  configuration <- as.Configuration(configuration)

  window <- as.owin(parameters$window)
  x_range <- window$xrange
  y_range <- window$yrange
  x_axis <- seq(from = x_range[1], to = x_range[2], length.out = grid_steps[1])
  y_axis <- seq(from = y_range[1], to = y_range[2], length.out = grid_steps[2])
  df <- as.data.frame(expand.grid(x = x_axis, y = y_axis))

  df$papangelou <- compute_papangelou_cpp(x = df$x,
                                          y = df$y,
                                          type = rep(type, length(df$x)),
                                          mark = rep(mark, length(df$x)),
                                          configuration = configuration,
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

  if(use_log) {
    df$papangelou <- log(df$papangelou)
    title <- "log-Papangelou conditional intensity"
  } else {
    title <- "Papangelou conditional intensity"
  }

  papangelou_im <-as.im(t(outer(x_axis, y_axis, function(x, y) {
    df$papangelou[df$x == x & df$y == y]
  })), W = window)

  if(return_papangelou) {
    return(papangelou_im)
  }

  if(use_ggplot) {
    points <- data.frame(x = configuration$x,
                         y = configuration$y,
                         types = droplevels(configuration$types),
                         marks = configuration$marks)
    color <- rep(c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"), nlevels(points$types))
    shape <- rep(c(16, 17, 15, 18), nlevels(points$types))

    if(missing(limits)) {
      lim <- c(min(df$papangelou), max(df$papangelou))
    } else {
      lim <- limits
    }

    g <- ggplot(data = df, aes_string(x = 'x', y = 'y')) +
      geom_tile(aes_string(fill = 'papangelou'), alpha = 0.5) +
      scale_fill_continuous_sequential(palette = "Purple-Yellow", limits = lim) +
      labs(fill = title) +
      scale_color_manual(values = color) + # Obtained by wesanderson::wes_palette("Darjeeling1")
      scale_shape_manual(values = shape) +
      xlab(NULL) + # Remove x labels
      ylab(NULL) + # Remove y labels
      ggtitle("") +
      coord_equal() +
      theme_minimal(base_size = 12) + # Theme
      guides(shape = guide_legend(override.aes = list(size = 5))) # Bigger species symbol size

    if(!all(points$marks == 1.)) {
      g <- g + geom_point(data = points, aes_string(x = 'x', y = 'y', colour = 'types', shape = 'types', size = 'marks'), alpha = 0.5)
    } else {
      g <- g + geom_point(data = points, aes_string(x = 'x', y = 'y', colour = 'types', shape = 'types'), size = 1.5, alpha = 0.5)
    }
    g
  } else {
    z <- outer(x_axis, y_axis, function(x, y) {
      df$papangelou[df$x == x & df$y == y]
    })
    plot(papangelou_im, main = title)
    plot(as.ppp(configuration, window), add = TRUE, cols = 'white')
  }
}

#' Plot the Papangelou conditional intensity from a fit object.
#'
#' @param fit Fit object obtained by running gibbsm.
#' @param window Observation window.
#' @param type Type of the point at which to evaluate the Papangelou conditional intensity.
#' @param configuration Configuration of points at which to evaluate the Papangelou conditional intensity.
#' @param alpha Short range repulsion parameter.
#' @param gamma Medium range repulsion parameter.
#' @param beta0 Log-intensity.
#' @param covariates Covariates.
#' @param beta Covariates coefficients.
#' @param short_range Matrix of short range distances.
#' @param medium_range Matrix of medium range distances.
#' @param long_range Matrix of long range distances.
#' @param saturation Saturation.
#' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
#' @param model Model for short range interaction.
#' @param medium_range_model Model for medium range interaction.
#' @param nthreads Maximum number of threads for parallel computing.
#' @param ... Forwarded to plot_papangelou
#' @export
#' @method plot_papangelou gibbsm
plot_papangelou.gibbsm <- function(fit,
                                   window = fit$window,
                                   type = 1,
                                   configuration = fit$configuration_list[[1]],
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
                                   nthreads = fit$nthreads,
                                   ...) {
  plot_papangelou(window = window,
                  type = type,
                  configuration = configuration,
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
                  nthreads = nthreads,
                  ...)
}
