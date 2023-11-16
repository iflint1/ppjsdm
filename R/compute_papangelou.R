#' Compute the Papangelou conditional intensity of the model.
#'
#' IMPORTANT: Check ?compute_papangelou.default for the documentation.
#'
#' @param ... Parameters to be forwarded to the relevant method. Currently, either
#' the parameters of the Gibbs point process, or a fit object obtained from running `gibbsm`.
#' @export
compute_papangelou <- function(...) {
  UseMethod("compute_papangelou")
}

#' Compute the Papangelou conditional intensity of the model with given parameters.
#'
#' @param configuration Configuration.
#' @param x Coordinates along the x-axis of the points at which to evaluate the Papangelou conditional intensity.
#' @param y Coordinates along the y-axis of the points at which to evaluate the Papangelou conditional intensity.
#' @param type Type of the point (as an integer >= 1 or string representing the type).
#' @param mark Mark of the point to add.
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
#' @param drop_type_from_configuration Should we remove the considered type(s) from the configuration?
#' @param nthreads Maximum number of threads for parallel computing.
#' @param ... Ignored.
#' @export
#' @method compute_papangelou default
compute_papangelou.default <- function(configuration,
                                       x,
                                       y,
                                       type = rep(1, length(x)),
                                       mark = rep(1.0, length(x)),
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
                                       drop_type_from_configuration = FALSE,
                                       nthreads = 1,
                                       ...) {
  # If user did not supply types, by default they should be those of the configuration
  if(missing(types) & !missing(configuration)) {
    if(nlevels(as.Configuration(configuration)$types) > 0) {
      types <- levels(as.Configuration(configuration)$types)
    } else {
      types <- NULL
    }

    # If the user supplied type, then add them to the list of types
    if(all(is.character(type))) {
      types <- sort(union(types, unique(type)))
    }
  }

  # Check format of type and mark parameters
  check_type_mark <- function(obj, force_numeric) {
    if(length(obj) != length(x)) {
      if(length(obj) == 1) {
        rep(obj, length(x))
      } else {
        stop("Unknown format for type or mark.")
      }
    } else {
      if(all(is.numeric(obj)) | (all(is.character(obj)) & !force_numeric)) {
        obj
      } else {
        stop("Type or mark are not of the right class.")
      }
    }
  }

  type <- check_type_mark(type, force_numeric = FALSE)
  mark <- check_type_mark(mark, force_numeric = TRUE)

  # Relevant types contains types mentioned by the user by name, we want to make sure that these are included
  relevant_types <- if(is.character(type)) {
    sort(unique(type))
  } else {
    NULL
  }

  if(!missing(configuration)) {
    configuration <- as.Configuration(configuration)
    if(nlevels(configuration$types) > 0) {
      relevant_types <- sort(union(relevant_types, levels(configuration$types)))
    }
  } else {
    configuration <- Configuration()
  }

  # Relevant indices contains types mentioned by the user by index
  relevant_indices <- if(is.numeric(type)) {
    unique(type)
  } else {
    NULL
  }

  # Construct defaults for the rest of the parameters
  parameters <- model_parameters(alpha = alpha,
                                 gamma = gamma,
                                 beta0 = beta0,
                                 covariates = covariates,
                                 beta = beta,
                                 short_range = short_range,
                                 medium_range = medium_range,
                                 long_range = long_range,
                                 saturation = saturation,
                                 model = model,
                                 types = types,
                                 medium_range_model = medium_range_model,
                                 relevant_types = relevant_types,
                                 relevant_indices = relevant_indices)

  # If the user supplies a string as the type, they want to evaluate the intensity
  # at that type.
  if(is.character(type)) {
    type <- sapply(type, function(t) which(t == names(parameters$beta0))[1])
  } else {
    type <- sapply(type, function(i) which(parameters$full_types[i] == names(parameters$beta0))[1])
  }

  # Remove type from the configuration?
  if(drop_type_from_configuration) {
    configuration <- configuration[setdiff(configuration$types, names(parameters$beta0)[unique(type)])]
  }

  # Check type format to avoid segfaults in C++ code
  if(!(is.integer(type) & all(type >= 1 & type <= length(parameters$beta0)))) {
    stop("Type does not have the right format")
  }

  # Fix type levels (configuration types might have fewer levels than the ones implied by other parameters,
  # causing issues in the C++ code)
  configuration$types <- factor(as.character(configuration$types), levels = names(parameters$beta0))

  compute_papangelou_cpp(x = x,
                         y = y,
                         type = type,
                         mark = mark,
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
}

#' Compute the Papangelou conditional intensity of the model from a fit object.
#'
#' @param fit Fit object obtained by running gibbsm.
#' @param type Type of the point at which to evaluate the Papangelou conditional intensity.
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
#' @param nthreads Maximum number of threads for parallel computing.
#' @param ... Forwarded to plot_papangelou
#' @export
#' @method compute_papangelou gibbsm
compute_papangelou.gibbsm <- function(fit,
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
  compute_papangelou(type = type,
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
