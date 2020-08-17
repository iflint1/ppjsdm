add_names <- function(str, covariates) {
  if(is.null(names(covariates))) {
    no_name <- rep(TRUE, length(covariates))
  } else {
    no_name <- names(covariates) == ""
  }
  names(covariates)[no_name] <- sprintf(paste0(str, "%d"), seq_len(length(which(no_name))))
  covariates
}

coerce_to_named_im_objects <- function(lst, str, window) {
  lst <- lapply(as.list(lst), function(element) {
    if(!is(element, "im")) {
      as.im(element, W = as.owin(window))
    } else {
      element
    }
  })
  add_names(str, lst)
}

model_parameters_defaults <- function(window,
                                      covariates,
                                      saturation,
                                      model,
                                      medium_range_model) {
  # The types below have default values.
  if(missing(window)) {
    window <- Rectangle_window()
  }
  if(missing(saturation)) {
    saturation <- 2
  }
  if(missing(model)) {
    model <- "square_bump"
  }
  if(missing(medium_range_model)) {
    medium_range_model <- "square_exponential"
  }
  if(missing(covariates)) {
    covariates <- list()
  } else {
    # Make covariates im objects with proper names.
    covariates <- coerce_to_named_im_objects(covariates, "unnamed_covariate", window)
  }

  list(window = window,
       saturation = saturation,
       model = model,
       medium_range_model = medium_range_model,
       covariates = covariates)
}

#' Construct parameters for the model.
#'
#' @param window Simulation window.
#' @param alpha Repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.
#' @param gamma Medium range repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.
#' @param beta0 A vector representing the log-intensities of the point processes.
#' Default is a vector of same size as types, filled with ones.
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
#' @param default_number_types Default number of types.
#' @export
model_parameters <- function(window,
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
                             default_number_types) {
  # The parameters below are set to NULL, and are treated in the C++ code.
  # Their values is not straightforward to determine and depends on how many types the user has in mind.
  if(missing(alpha)) {
    alpha <- NULL
  }

  if(missing(gamma)) {
    gamma <- NULL
  }

  if(missing(beta0)) {
    beta0 <- NULL
  }

  if(missing(beta)) {
    beta <- NULL
  } else {
    beta <- as.matrix(beta)
  }

  if(missing(short_range)) {
    short_range <- NULL
  }

  if(missing(medium_range)) {
    medium_range <- NULL
  }

  if(missing(long_range)) {
    long_range <- NULL
  }

  if(missing(types)) {
    types <- NULL
  }

  if(missing(default_number_types)) {
    default_number_types <- 1
  }

  defaults <- model_parameters_defaults(window = window,
                                        covariates = covariates,
                                        saturation = saturation,
                                        model = model,
                                        medium_range_model = medium_range_model)

  parameters <- make_default_model_parameters(alpha = alpha,
                                              beta0 = beta0,
                                              covariates = defaults$covariates,
                                              beta = beta,
                                              gamma = gamma,
                                              short_range = short_range,
                                              medium_range = medium_range,
                                              long_range = long_range,
                                              types = types,
                                              default_number_types = default_number_types)
  append(parameters, defaults)
}
