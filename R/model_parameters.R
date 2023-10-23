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
  # Model is a list of length the number of potentials, missing values indicated by NULL
  model <- lapply(model, function(x) {
    if(is.null(x)) {
      "square_bump"
    } else {
      x
    }
  })
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

#' Construct an object that might already exist as `x`, but if not that should be constructed
#' using default supplied values, as either a matrix or a vector.
#'
#' @param x Object that contains either an object, or NULL if it should be constructed.
#' @param def Default value for the object.
#' @param nrows Number of rows of the object to construct.
#' @param ncols Number of columns of the object to construct.
#' @param matrix Should we construct a matrix or a vector?
#' @export
#' @keywords internal
construct_if_missing <- function(x, def, nrows, ncols = nrows, matrix = FALSE) {
  if(is.null(x)) {
    if(matrix) {
      matrix(def, ncol = ncols, nrow = nrows)
    } else {
      rep(def, times = nrows)
    }
  } else if(matrix) {
    if(!is.matrix(x)) {
      if(is.numeric(x) & length(x) == 1) {
        matrix(x, ncol = ncols, nrow = nrows)
      } else {
        stop(paste0("Not sure how to convert one of the parameters to a matrix: ", x, sep = ", "))
      }
    } else {
      as.matrix(x)
    }
  } else {
    setNames(as.numeric(x), nm = names(x))
  }
}

get_number_types_implied_by_object <- function(x) {
  if(is.null(x)) {
    0
  } else if(is.matrix(x)) {
    length_x <- nrow(x)
    if(length_x != ncol(x)) {
      stop("Found a non-square matrix in the arguments.")
    }
    length_x
  } else if(is.numeric(x) | is.character(x)) {
    length(x)
  } else if(is.list(x)) {
    initial_length <- get_number_types_implied_by_object(x[[1]])
    if(length(x) > 1) {
      for(i in 2:length(x)) {
        if(get_number_types_implied_by_object(x[[i]]) != initial_length) {
          stop("At least two objects in a list had incompatible sizes.")
        }
      }
    }
    initial_length
  } else {
    stop(paste0("Not sure how to get number of types from this object: ", x))
  }
}

get_number_types_and_check_conformance_helper <- function(default_number_types,
                                                          best_guess,
                                                          x,
                                                          ...) {
  if(missing(x)) {
    if(best_guess == 0) {
      default_number_types
    } else {
      best_guess
    }
  } else {
    length_x <- get_number_types_implied_by_object(x)
    if(is.null(x) | length_x == 0) {
      get_number_types_and_check_conformance_helper(default_number_types = default_number_types,
                                                    best_guess = best_guess,
                                                    ...)
    } else {
      if(best_guess != 0) {
        if(length_x != best_guess && (best_guess > 1 && length_x > 1)) {
          stop("Two of the given arguments have incompatible sizes.")
        } else {
          get_number_types_and_check_conformance_helper(default_number_types = default_number_types,
                                                        best_guess = max(best_guess, length_x),
                                                        ...)
        }
      } else {
        get_number_types_and_check_conformance_helper(default_number_types = default_number_types,
                                                      best_guess = length_x,
                                                      ...)
      }
    }
  }
}

#' Deduce from a list of objects the number of types the user has in mind.
#'
#' @param default_number_types Default number of types.
#' @param ... Objects (usually matrices and vectors) that imply varying number of types.
#' @export
#' @keywords internal
get_number_types_and_check_conformance <- function(default_number_types,
                                                   ...) {
  get_number_types_and_check_conformance_helper(default_number_types = default_number_types,
                                                best_guess = 0,
                                                ...)
}

make_default_model_parameters <- function(alpha,
                                          beta0,
                                          covariates,
                                          beta,
                                          gamma,
                                          short_range,
                                          medium_range,
                                          long_range,
                                          types,
                                          default_number_types,
                                          ...) {
  number_types <- get_number_types_and_check_conformance(default_number_types = default_number_types,
                                                         alpha,
                                                         gamma,
                                                         beta0,
                                                         short_range,
                                                         medium_range,
                                                         long_range,
                                                         types,
                                                         ...)

  alpha <- lapply(alpha, function(a) {
    a <- construct_if_missing(a, 0, number_types, matrix = TRUE)
    if(!isSymmetric(a)) {
      stop("Alpha is not symmetric.");
    }
    a
  })

  gamma <- construct_if_missing(gamma, 0, number_types, matrix = TRUE)
  if(!isSymmetric(gamma)) {
    stop("Gamma is not symmetric.");
  }

  beta0 <- construct_if_missing(beta0, 0, number_types, matrix = FALSE)

  beta_nrows <- number_types
  beta_ncols <- length(covariates)

  beta = construct_if_missing(beta, 1., beta_nrows, beta_ncols, matrix = TRUE)
  if(beta_ncols != 0 && (nrow(beta) != beta_nrows || ncol(beta) != beta_ncols)) {
    stop("The parameter `beta` does not have the right dimensions.");
  }

  types <- make_types(types = types,
                      size = number_types,
                      might_contain_name = beta0)

  default_short_distances <- seq(from = 0, to = 0.1, length.out = length(short_range) + 1)[-1]
  short_range <- lapply(seq_len(length(short_range)), function(i) {
    s <- construct_if_missing(short_range[[i]], default_short_distances[i], number_types, matrix = TRUE)
    if(!isSymmetric(s)) {
      stop("One of the short-range interaction radii matrices is not symmetric.");
    }
    s
  })

  medium_range <- construct_if_missing(medium_range, 0, number_types, matrix = TRUE)
  if(!isSymmetric(medium_range)) {
    stop("The medium-range interaction radii matrix is not symmetric.");
  }

  long_range <- construct_if_missing(long_range, 0, number_types, matrix = TRUE)
  if(!isSymmetric(long_range)) {
    stop("The long-range interaction radii matrix is not symmetric.");
  }

  list(alpha = alpha,
       beta0 = beta0,
       covariates = covariates,
       beta = beta,
       gamma = gamma,
       short_range = short_range,
       medium_range = medium_range,
       long_range = long_range,
       types = types,
       ...);
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
#' @param ... Other parameters used to infer the number of types.
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
                             default_number_types = 1,
                             ...) {
  # TODO: A lot of this is converted from C++ code, and has been tested by hand on a few scenarios.
  # WRITE PROPER TESTS to make sure it all works as expected. E.g., numeric short_range, varying number of potentials specified
  # via short_range, model, alpha, etc.

  # Try to figure out how many short_range potentials the user intends to use
  number_short_range <- 1
  if(!missing(alpha)) {
    if(is.list(alpha)) {
      number_short_range <- length(alpha)
    }
  }
  if(!missing(short_range)) {
    if(is.list(short_range)) {
      if(number_short_range > 1 & number_short_range != length(short_range)) {
        stop("Alpha was a list, indicating >1 potentials, but short_range does not have the same number of potentials.")
      }
      number_short_range <- length(short_range)
    }
  }
  if(!missing(model)) {
    if(is.list(model)) {
      if(number_short_range > 1 & number_short_range != length(model)) {
        stop("Either alpha or short_range was a list, indicating >1 potentials, but model does not have the same number of potentials.")
      }
      number_short_range <- length(model)
    }
  }

  # Default initialise most parameters, but at this point we do not know how many types user intends
  if(missing(alpha)) {
    alpha <- lapply(seq_len(number_short_range), function(x) NULL)
  } else if(!is.list(alpha)) {
    if(is.numeric(alpha)) {
      alpha <- lapply(seq_len(number_short_range), function(x) alpha)
    } else {
      stop("Not sure how to interpret parameter alpha.")
    }
  }

  if(missing(short_range)) {
    short_range <- lapply(seq_len(number_short_range), function(x) NULL)
  } else if(!is.list(short_range)) {
    if(is.numeric(short_range)) {
      if(number_short_range > 1) {
        stop("Detected a numeric short_range (implying a single short-range potential) but other parameters indicate >1 short-range potentials. There is likely a mistake in the model specification.")
      }
      short_range <- lapply(seq_len(number_short_range), function(x) short_range)
    } else {
      stop("Not sure how to interpret parameter short_range")
    }
  }

  if(missing(model)) {
    model <- lapply(seq_len(number_short_range), function(x) NULL)
  } else if(!is.list(short_range)) {
    if(is.character(short_range)) {
      model <- lapply(seq_len(number_short_range), function(x) model)
    } else {
      stop("Not sure how to interpret parameter model")
    }
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

  if(missing(medium_range)) {
    medium_range <- NULL
  }

  if(missing(long_range)) {
    long_range <- NULL
  }

  if(missing(types)) {
    types <- NULL
  } else {
    types <- unlist(types)
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
                                              default_number_types = default_number_types,
                                              ...)

  if(!is.list(defaults$model)) {
    defaults$model <- lapply(parameters$alpha, function(x) defaults$model)
  }

  if(length(defaults$model) != length(parameters$alpha) || length(defaults$model) != length(parameters$short_range)) {
    stop(paste0("Some of the parameters have incompatible sizes; model: [", paste0(defaults$model, collapse = ", "),
                "], alpha: [", paste0(parameters$alpha, collapse = ", "),
                "] and short_range: [", paste0(parameters$short_range, collapse = ", "), "]."))
  }

  parameters$alpha <- lapply(parameters$alpha, function(a) {
    z <- a
    colnames(z) <- rownames(z) <- parameters$types
    z
  })
  colnames(parameters$gamma) <- rownames(parameters$gamma) <- parameters$types
  parameters$short_range <- lapply(parameters$short_range, function(s) {
    z <- s
    colnames(z) <- rownames(z) <- parameters$types
    z
  })
  colnames(parameters$medium_range) <- rownames(parameters$medium_range) <- parameters$types
  colnames(parameters$long_range) <- rownames(parameters$long_range) <- parameters$types
  rownames(parameters$beta) <- parameters$types
  names(parameters$beta0) <- parameters$types
  append(parameters, defaults)
}
