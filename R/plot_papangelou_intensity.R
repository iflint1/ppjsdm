#' Plot the Papangelou conditional intensity.
#'
#' IMPORTANT: Check ?plot_papangelou.default for the documentation. Here, we only provide a few examples.
#'
#' @param ... Parameters to be forwarded to the relevant method. Currently, either
#' the parameters of the Gibbs point process, or a fit object obtained from running `gibbsm`.
#' @export
#' @examples
#' # Plot the Papangelou conditional intensity of a model specified by some parameters.
#' # This will generate a configuration from that model, and plot the corresponding intensity.
#'
#' ppjsdm::plot_papangelou(alpha = matrix(-0.2, 2, 2), # Short-range interaction coefficients
#'                         beta0 = c(4, 4),
#'                         grid_steps = 100)
#'
#' # More commonly, one starts from a set of points:
#'
#' configuration <- ppjsdm::rppp(lambda = c(A = 500, B = 500))
#'
#' # Which is used to fit the model:
#'
#' fit <- ppjsdm::gibbsm(configuration)
#'
#' # And the fit is then used to generate the conditional intensity:
#'
#' ppjsdm::plot_papangelou(fit, type = "A", grid_steps = 100)
#'
#' # It is possible to override any of the model parameters, e.g., the window:
#'
#' ppjsdm::plot_papangelou(fit, type = "A", grid_steps = 100, window = ppjsdm::Rectangle_window(c(0, 0.5), c(0, 1)))
#'
#' @md
plot_papangelou <- function(...) {
  UseMethod("plot_papangelou")
}

#' Plot the Papangelou conditional intensity with given parameters.
#'
#' @param window Observation window. Preferably a `ppjsdm` Window, such as `ppjsdm::Rectangle_window`, but also accepts `spatstat` `im` or `owin` objects.
#' @param type Type of the point at which to evaluate the Papangelou conditional intensity. Either an integer representing the index of the type, or a string representing its label.
#' @param mark Mark of the point.
#' @param configuration Configuration of points at which to evaluate the Papangelou conditional intensity.
#' @param alpha Matrix of short-range interaction coefficients.
#' @param gamma Matrix of medium-range interaction coefficients.
#' @param beta0 Vector representing the intercept.
#' @param covariates List of covariates. These are converted to the `im` format by applying `as.im` to all elements in the list.
#' @param beta Fitted regression coefficients with respect to covariates.
#' @param short_range Symmetric matrix of short-range interaction radii. Can also be a list of matrices, each entry representing a different potential.
#' @param medium_range Symmetric matrix of medium-range interaction radii.
#' @param long_range Symmetric matrix of long-range interaction radii.
#' @param saturation Saturation parameter of the point process.
#' @param types Character vector, with entry i representing the name of type i.
#' @param model String representing the short-range model to use. The currently authorised models are obtained with a call to `show_short_range_models()`.
#' @param medium_range_model String representing the medium-range model to use. The currently authorised models are obtained with a call to `show_medium_range_models()`.
#' @param grid_steps Vector of length 2, representing the number of horizontal and vertical grid steps. If a single value is given, it is assumed that the dimensions are the same along both axes.
#' @param steps Number of steps in the Metropolis-Hastings simulation algorithm. A value of 0 uses the CFTP algorithm.
#' This is only used if no configuration was supplied, in which case we simulate a configuration.
#' @param nthreads Maximum number of threads to use.
#' @param use_log Plot the logarithm of the Papangelou conditional intensity instead?
#' @param use_ggplot Use ggplot for fancier plot?
#' @param return_papangelou Should we return the Papangelou intensity as an `im` object instead of plotting it?
#' @param drop_type_from_configuration Should we remove the considered `type` from the configuration?
#' @param show In the plot, should we show all points, no points, or only the individuals of the configuration with the target type?
#' @param limits Limits for values of the Papangelou conditional intensity when plotting.
#' @param ... Ignored.
#' @importFrom colorspace scale_fill_continuous_sequential
#' @importFrom ggplot2 aes coord_equal element_text geom_tile ggplot ggtitle guide_legend guides labs scale_color_manual scale_shape_manual scale_x_continuous scale_y_continuous theme theme_minimal xlab xlim ylab ylim
#' @importFrom graphics plot
#' @importFrom spatstat.geom as.im as.owin as.ppp boundingbox gridcentres inside.owin
#' @importFrom stats na.omit
#' @export
#' @method plot_papangelou default
#' @md
plot_papangelou.default <- function(window,
                                    type = 1,
                                    mark = 1.0,
                                    configuration,
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
                                    grid_steps = c(1000, 1000),
                                    steps = 0,
                                    nthreads = 1,
                                    use_log = FALSE,
                                    use_ggplot = TRUE,
                                    return_papangelou = FALSE,
                                    drop_type_from_configuration = FALSE,
                                    show = c("all", "none", "type"),
                                    limits,
                                    ...) {
  # Interpret the show argument
  show <- match.arg(show)

  # Take care of grid_steps argument
  if(!is.numeric(grid_steps)) {
    stop("Expected numeric grid_steps")
  } else {
    if(length(grid_steps) == 1) {
      grid_steps <- c(grid_steps, grid_steps)
    } else if(length(grid_steps) != 2) {
      stop("Expected grid_steps to have length 1 or 2.")
    }
  }

  # If user did not supply types, by default they should be those of the configuration
  if(missing(types) & !missing(configuration)) {
    types <- levels(as.Configuration(configuration)$types)
  }

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
                            steps = steps,
                            nthreads = nthreads)
  }
  configuration <- as.Configuration(configuration)

  # TODO: There is a lot of code duplication here and in compute_papangelou, try to factorise...
  # There are some cases where configuration might only have a subset of all types,
  # e.g., if one wants to predict at a new location.
  if(nlevels(configuration$types) != length(parameters$types)) {
    # In this case, everything is ok, we can evaluate the parameters on the subset of types
    if(all(levels(configuration$types) %in% parameters$types)) {
      parameters$alpha <- lapply(parameters$alpha, function(a) a[levels(configuration$types), levels(configuration$types)])
      parameters$gamma <- parameters$gamma[levels(configuration$types), levels(configuration$types)]
      parameters$short_range <- lapply(parameters$short_range, function(s) s[levels(configuration$types), levels(configuration$types)])
      parameters$medium_range <- parameters$medium_range[levels(configuration$types), levels(configuration$types)]
      parameters$long_range <- parameters$long_range[levels(configuration$types), levels(configuration$types)]

      parameters$beta0 <- parameters$beta0[levels(configuration$types)]
      parameters$beta <- parameters$beta[levels(configuration$types), ]

      parameters$types <- levels(configuration$types)
    } else {
      stop(paste0("The types of the configuration are not a subset of those given by the parameters, configuration: ",
                  paste0(levels(configuration$types), collapse = ", "),
                  " and supplied types: ",
                  paste0(parameters$types, collapse = ", ")))
    }
  } else if(!all(levels(configuration$types) == parameters$types)) {
    stop(paste0("The types of the configuration do not correspond to those given by the parameters, configuration: ",
                paste0(levels(configuration$types), collapse = ", "),
                " and supplied types: ",
                paste0(parameters$types, collapse = ", ")))
  }

  # At this point, the parameters should exactly correspond to the types of the configuration
  if(!all(names(parameters$beta0) == levels(configuration$types))) {
    stop(paste0("Unexpected setting when computing the Papangelou intensity: the supplied parameters do not refer to the same types as the configuration, configuration: ",
                paste0(levels(configuration$types), collapse = ", "),
                " and supplied parameters: ",
                paste0(names(parameters$beta0), collapse = ", ")))
  }

  # If the user supplies a string as the type, they want to evaluate the intensity
  # at that type.
  if(length(type) != 1) {
    stop("Expecting a type of length 1.")
  } else if(is.character(type)) {
    type <- which(type == names(parameters$beta0))[1]
  } else if(is.numeric(type)) {
    type <- as.integer(type)
  } else {
    stop("Expecting a type that is either an integer or a character.")
  }

  # Take care of mark and type arguments
  if(!is.numeric(mark) | length(mark) != 1) {
    stop("Expecting mark to be a single numeric value.")
  }

  window <- as.owin(parameters$window)
  df <- as.data.frame(gridcentres(window, nx = grid_steps[1], ny = grid_steps[2]))
  df <- df[inside.owin(x = df$x, y = df$y, w = window), ]

  df$papangelou <- compute_papangelou(x = df$x,
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
                                      types = parameters$types,
                                      drop_type_from_configuration = drop_type_from_configuration,
                                      nthreads = nthreads)

  if(use_log) {
    df$papangelou <- log(df$papangelou)
    title <- "log-Papangelou conditional intensity"
  } else {
    title <- "Papangelou conditional intensity"
  }

  papangelou_im <- as.im(df)

  if(return_papangelou) {
    return(papangelou_im)
  }

  # Which individuals do we want to show on the plot?
  if(show == "type") {
    configuration <- Configuration(x = configuration$x[configuration$types %in% parameters$types[type]],
                                   y = configuration$y[configuration$types %in% parameters$types[type]],
                                   types = configuration$types[configuration$types %in% parameters$types[type]],
                                   marks = configuration$marks[configuration$types %in% parameters$types[type]])
  } else if(show == "none") {
    configuration <- Configuration()
  }

  if(use_ggplot) {
    points <- data.frame(x = configuration$x,
                         y = configuration$y,
                         Types = droplevels(configuration$types),
                         marks = configuration$marks)
    points <- points[inside.owin(x = points$x, y = points$y, w = window), ]
    color <- rep(c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"), nlevels(points$Types))
    shape <- rep(c(16, 17, 15, 18), nlevels(points$Types))

    lim <- if(missing(limits)) {
      c(min(df$papangelou), max(df$papangelou))
    } else {
      limits
    }

    g <- ggplot(data = df, aes_string(x = "x", y = "y")) +
      geom_tile(aes_string(fill = "papangelou"), alpha = 0.5) +
      scale_fill_continuous_sequential(palette = "Purple-Yellow", limits = lim) +
      labs(fill = title) +
      scale_color_manual(values = color) + # Obtained by wesanderson::wes_palette("Darjeeling1")
      scale_shape_manual(values = shape) +
      xlab(NULL) + # Remove x labels
      ylab(NULL) + # Remove y labels
      ggtitle("") +
      coord_equal() +
      theme_minimal(base_size = 12) + # Theme
      xlim(boundingbox(window)$xrange) +
      ylim(boundingbox(window)$yrange) +
      guides(shape = guide_legend(override.aes = list(size = 5))) # Bigger species symbol size

    if(!all(points$marks == 1.)) {
      g <- g + geom_point(data = points, aes_string(x = 'x', y = 'y', colour = 'Types', shape = 'Types', size = 'marks'), alpha = 0.5)
    } else {
      g <- g + geom_point(data = points, aes_string(x = 'x', y = 'y', colour = 'Types', shape = 'Types'), size = 1.5, alpha = 0.5)
    }
    g
  } else {
    plot(papangelou_im, main = title)
    plot(as.ppp(configuration, window), add = TRUE, cols = 'white')
  }
}

#' Plot the Papangelou conditional intensity from a fit object.
#'
#' @param fit Fit object obtained by running gibbsm.
#' @param window Observation window. Preferably a `ppjsdm` Window, such as `ppjsdm::Rectangle_window`, but also accepts `spatstat` `im` or `owin` objects.
#' @param configuration Configuration of points at which to evaluate the Papangelou conditional intensity.
#' @param type Type of the point at which to evaluate the Papangelou conditional intensity. Either an integer representing the index of the type, or a string representing its label.
#' @param alpha Matrix of short-range interaction coefficients.
#' @param gamma Matrix of medium-range interaction coefficients.
#' @param beta0 Vector representing the intercept.
#' @param covariates List of covariates. These are converted to the `im` format by applying `as.im` to all elements in the list.
#' @param beta Fitted regression coefficients with respect to covariates.
#' @param short_range Symmetric matrix of short-range interaction radii. Can also be a list of matrices, each entry representing a different potential.
#' @param medium_range Symmetric matrix of medium-range interaction radii.
#' @param long_range Symmetric matrix of long-range interaction radii.
#' @param saturation Saturation parameter of the point process.
#' @param types Character vector, with entry i representing the name of type i.
#' @param model String representing the short-range model to use. The currently authorised models are obtained with a call to `show_short_range_models()`.
#' @param medium_range_model String representing the medium-range model to use. The currently authorised models are obtained with a call to `show_medium_range_models()`.
#' @param nthreads Maximum number of threads to use.
#' @param ... Forwarded to plot_papangelou
#' @export
#' @method plot_papangelou gibbsm
#' @md
plot_papangelou.gibbsm <- function(fit,
                                   window = fit$window,
                                   configuration = fit$configuration_list[[1]],
                                   type = 1,
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
                  configuration = configuration,
                  type = type, # Need to have this here and NOT forwarded, otherwise a given type argument gets matched with types
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
