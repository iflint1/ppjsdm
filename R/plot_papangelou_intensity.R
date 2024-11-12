#' Plot the Papangelou conditional intensity.
#'
#' IMPORTANT: Check ?plot_papangelou.default for the documentation. Here, we only provide a few examples.
#'
#' @param ... Parameters to be forwarded to the relevant method. Currently, either
#' the parameters of the Gibbs point process, or a fit object obtained from running `gibbsm`.
#' @export
#' @examples
#' set.seed(1)
#'
#' # Draw some points
#'
#' configuration <- ppjsdm::rppp(lambda = c(A = 500, B = 500))
#'
#' # Plot the Papangelou conditional intensity of a model specified by some parameters.
#'
#' ppjsdm::plot_papangelou(configuration = configuration,
#'                         alpha = matrix(-0.2, 2, 2), # Short-range interaction coefficients
#'                         beta0 = c(4, 4),
#'                         grid_steps = 100)
#'
#' # More commonly, the points are used to fit the model:
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
#' @param nthreads Maximum number of threads to use.
#' @param use_log Plot the logarithm of the Papangelou conditional intensity instead?
#' @param use_ggplot Use ggplot for fancier plot?
#' @param return_papangelou Should we return the Papangelou intensity as an `im` object instead of plotting it?
#' @param drop_type_from_configuration Should we remove the considered `type` from the configuration?
#' @param show In the plot, should we show all points, no points, or only the individuals of the configuration with certain types?
#' `show = "all"` shows all points, `show = "none"` shows no points, `show = "type"` shows only the focal type,
#' and `show` can otherwise be a vector of types, in which case all the corresponding types will be shown.
#' @param limits Limits for values of the Papangelou conditional intensity when plotting.
#' @param mark_range (Optional) Vector of length two, representing the range of sizes used when plotting points with different marks.
#' @param base_size Font size.
#' @param full_configuration (Optional) Configuration object, perhaps different from `configuration`, to overlay on the Papangelou conditional intensity.
#' This need not be the configuration used to compute the intensity, indeed, one might wanr to condition on some types, but show others on the final plot.
#' @param legend_title Title to give to the legend describing the conditional intensity.
#' @param type_description Name to give the types in the legend (e.g., "species" or "classes").
#' @param mark_description Name to give the marks in the legend (e.g., "size" or "DBH").
#' @param types_order (Optional) Order of the types in figure legend. Vector of the names of types to plot, ordered in a given way.
#' @param colours (Optional) Colours to use when plotting points.
#' @param shapes (Optional) Shapes to use when plotting points.
#' @param ... Ignored.
#' @importFrom colorspace scale_fill_continuous_sequential
#' @importFrom ggplot2 aes coord_equal element_text geom_tile ggplot ggtitle guide_colourbar guide_legend guides labs scale_color_manual scale_shape_manual scale_size scale_x_continuous scale_y_continuous theme theme_minimal xlab xlim ylab ylim
#' @importFrom graphics plot
#' @importFrom rlang .data
#' @importFrom scales breaks_extended
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
                                    grid_steps = c(100, 100),
                                    nthreads = 1,
                                    use_log = FALSE,
                                    use_ggplot = TRUE,
                                    return_papangelou = FALSE,
                                    drop_type_from_configuration = FALSE,
                                    show = "all",
                                    limits,
                                    mark_range = c(1, 6),
                                    base_size = 11,
                                    full_configuration,
                                    legend_title,
                                    type_description = "Types",
                                    mark_description = "Marks",
                                    types_order,
                                    colours,
                                    shapes,
                                    ...) {

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

  # Take care of mark and type arguments
  if(!is.numeric(mark) | length(mark) != 1) {
    stop("Expecting mark to be a single numeric value.")
  }

  if(!(is.numeric(type) | is.character(type)) | length(type) != 1) {
    stop("Expecting a char/numeric type of length 1.")
  }

  # Convert window to owin, and subset grid points
  window <- as.owin(model_parameters(window = window)$window)
  df <- as.data.frame(gridcentres(window, nx = grid_steps[1], ny = grid_steps[2]))
  df <- df[inside.owin(x = df$x, y = df$y, w = window), ]

  # Subset configuration to window
  if(missing(configuration)) {
    configuration <- Configuration()
  } else {
    configuration <- as.data.frame(as.Configuration(configuration))
    configuration <- configuration[inside.owin(x = configuration$x,
                                               y = configuration$y,
                                               w = window), ]
    configuration <- as.Configuration(configuration)
  }

  # Subset full_configuration to window
  if(missing(full_configuration)) {
    full_configuration <- configuration
  } else {
    full_configuration <- as.data.frame(as.Configuration(full_configuration))
    full_configuration <- full_configuration[inside.owin(x = full_configuration$x,
                                                         y = full_configuration$y,
                                                         w = window), ]
    full_configuration <- as.Configuration(full_configuration)
  }

  df$papangelou <- compute_papangelou(x = df$x,
                                      y = df$y,
                                      type = type,
                                      mark = mark,
                                      configuration = configuration,
                                      model = model,
                                      medium_range_model = medium_range_model,
                                      alpha = alpha,
                                      beta0 = beta0,
                                      beta = beta,
                                      gamma = gamma,
                                      covariates = covariates,
                                      short_range = short_range,
                                      medium_range = medium_range,
                                      long_range = long_range,
                                      saturation = saturation,
                                      types = types,
                                      drop_type_from_configuration = drop_type_from_configuration,
                                      nthreads = nthreads)

  if(use_log) {
    df$papangelou <- log(df$papangelou)
    if(missing(legend_title)) {
      legend_title <- "log-Papangelou conditional intensity"
    }
  } else if(missing(legend_title)) {
    legend_title <- "Papangelou conditional intensity"
  }

  papangelou_im <- as.im(df)

  if(return_papangelou) {
    return(papangelou_im)
  }

  # Which individuals do we want to show on the plot?
  if(all(show %in% levels(full_configuration$types))) {
    full_configuration <- Configuration(x = full_configuration$x[full_configuration$types %in% show],
                                        y = full_configuration$y[full_configuration$types %in% show],
                                        types = full_configuration$types[full_configuration$types %in% show],
                                        marks = full_configuration$marks[full_configuration$types %in% show])
  } else if(length(show) == 1) {
    if(show == "none") {
      full_configuration <- Configuration()
    } else if(show == "type") {
      #full_configuration <- full_configuration[names(parameters$beta0)[type]]
      stop("Parameter show == \"type\" not allowed for now, use the actual name of the type to show instead.")
    } else if(show != "all") {
      print(show)
      stop("Could not interpret the parameter show, printed above.")
    }
  } else {
    print(show)
    stop("Parameter show was of length != 1, but its values were not types found in the configuration. The provided parameters was printed above.")
  }

  if(use_ggplot) {
    points <- as.data.frame(full_configuration)
    names(points)[names(points) == "types"] <- type_description
    names(points)[names(points) == "marks"] <- mark_description

    lim <- if(missing(limits)) {
      c(min(df$papangelou), max(df$papangelou))
    } else if(!is.null(limits) & !(is.numeric(limits) & length(limits) == 2)) {
      print(limits)
      stop("Did not understand supplied limits parameter (printed above)")
    } else {
      limits
    }

    if(missing(shapes)) {
      shapes <- rep(c(16, 17, 15, 18), length.out = nlevels(points[, type_description]))
    }

    if(!missing(types_order)) {
      points[, type_description] <- droplevels(points[, type_description])
      if(!identical(sort(types_order), sort(levels(points[, type_description])))) {
        stop("The parameter \"types_order\" should contain exactly the same name as those to be plotted, with perhaps their order changed.")
      }
      points[, type_description] <- factor(points[, type_description], levels = types_order)
    }

    g <- ggplot(data = df, aes(x = .data$x, y = .data$y)) +
      geom_tile(aes(fill = .data$papangelou), alpha = 1) +
      scale_fill_continuous_sequential(palette = "Purple-Yellow", limits = lim) +
      labs(fill = legend_title) +
      scale_shape_manual(values = shapes) +
      xlab(NULL) + # Remove x labels
      ylab(NULL) + # Remove y labels
      ggtitle("") +
      coord_equal() +
      theme_minimal(base_size = base_size) + # Theme
      xlim(boundingbox(window)$xrange) +
      ylim(boundingbox(window)$yrange)

    if(!all(points[, mark_description] == 1.)) {
      g <- g + geom_point(data = points, aes(x = .data$x, y = .data$y, colour = .data[[type_description]],
                                             shape = .data[[type_description]],
                                             size = .data[[mark_description]]), alpha = 0.8) +
        scale_size(name = mark_description, range = mark_range, breaks = breaks_extended(6))

      nr <- 6
    } else {
      g <- g + geom_point(data = points, aes(x = .data$x, y = .data$y, colour = .data[[type_description]],
                                             shape = .data[[type_description]]), size = mean(mark_range), alpha = 0.8)

      nr <- 8
    }

    if(!missing(colours)) {
      g <- g + scale_colour_manual(values = rep(colours, nlevels(points[, type_description])))
    } else {
      g <- g + scale_colour_viridis_d(end = 0.9, option = "turbo")
    }

    if(!all(points[, mark_description] == 1.)) {
      g <- g + guides(colour = guide_legend(order = 1,
                                            nrow = nr,
                                            override.aes = list(size = 5)),
                      shape = guide_legend(order = 1,
                                           nrow = nr),
                      size = guide_legend(order = 2),
                      fill = guide_colourbar(order = 3))
    } else {
      g <- g + guides(colour = guide_legend(order = 1,
                                            nrow = nr,
                                            override.aes = list(size = 5)),
                      shape = guide_legend(order = 1,
                                           nrow = nr),
                      fill = guide_colourbar(order = 2))
    }

    g
  } else {
    plot(papangelou_im, main = legend_title)
    plot(as.ppp(full_configuration, window), add = TRUE, cols = 'white')
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
#' @param full_configuration (Optional) Configuration object, perhaps different from `configuration`, to overlay on the Papangelou conditional intensity.
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
                                   full_configuration = fit$configuration_list[[1]],
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
                  full_configuration = full_configuration,
                  ...)
}
