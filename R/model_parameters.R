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

#' Construct an object `x` that might already be constructed, but if not that should be constructed
#' using default supplied values, as either a matrix or a vector.
#'
#' @param x Object that contains either an object, or NULL if it should be constructed.
#' @param def Default value for the object.
#' @param nrows Number of rows of the object to construct.
#' @param ncols Number of columns of the object to construct.
#' @param matrix Should we construct a matrix or a vector?
#' @export
#' @keywords internal
#' @md
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

get_number_types_implied_by_object <- function(x, types) {
  # TODO: Horrible function, factorise...
  if(is.null(x)) {
    list(length = 0,
         types = types)
  } else if(is.matrix(x)) {
    if(is.null(types)) {
      ty <- if(!is.null(rownames(x))) {
        rownames(x)
      } else {
        NULL
      }
      list(length = nrow(x),
           types = ty)
    } else if(!is.null(rownames(x))) {
      if(all(types %in% rownames(x))) { # The matrix x will be subset to only include the relevant types
        list(length = nrow(x),
             types = rownames(x))
      } else if(all(rownames(x) %in% types)) { # The types will be subset to only include the rows in x
        list(length = length(types),
             types = types)
      } else {
        print(x)
        print(types)
        stop("Found a matrix (printed above) that has rownames, but provided types (printed afterwards) that are not compatible.")
      }
    } else {
      list(length = nrow(x),
           types = types)
    }
  } else if(is.numeric(x) | is.character(x)) {
    if(is.null(types)) {
      ty <- if(!is.null(names(x))) {
        names(x)
      } else {
        NULL
      }
      list(length = length(x),
           types = ty)
    } else if(!is.null(names(x))) {
      if(all(types %in% names(x))) { # The vector x will be subset to only include the relevant types
        list(length = length(x),
             types = names(x))
      } else if(all(names(x) %in% types)) {
        list(length = length(types),
             types = types)
      } else {
        print(x)
        print(types)
        stop("Found a vector (printed above) that has rownames, but provided types (printed afterwards) that are not compatible.")
      }
    } else {
      list(length = length(x),
           types = types)
    }
  } else if(is.list(x)) {
    g <- get_number_types_implied_by_object(x[[1]], types = types)
    initial_length <- g$length
    ty <- g$types
    if(!all(types %in% ty)) {
      stop(paste0("Considering an object ", x[[1]], " that has implied type names incompatible with ", types))
    }
    if(length(x) > 1) {
      for(i in 2:length(x)) {
        g <- get_number_types_implied_by_object(x[[i]], types = ty)
        if(initial_length > 1 & g$length > 1 & g$length != initial_length) {
          stop("At least two objects in a list had incompatible sizes.")
        }
        initial_length <- max(initial_length, g$length)
        if(!all(ty %in% g$ty)) {
          stop(paste0("Considering an object ", x[[i]], " that has implied type names incompatible with ", ty))
        }
        ty <- g$ty
      }
    }
    list(length = initial_length,
         types = ty)
  } else {
    stop(paste0("Not sure how to get number of types from this object: ", x))
  }
}

#' Deduce from a list of objects the number of types the user has in mind.
#'
#' @param default_number_types Default number of types.
#' @param ... Objects (usually matrices and vectors) that imply varying number of types.
#' @export
#' @keywords internal
get_number_types_and_check_conformance <- function(default_number_types,
                                                   types,
                                                   ...) {
  best_guess <- if(!is.null(types)) {
    length(types)
  } else {
    0
  }

  best_types <- types

  for(x in list(...)) {
    g <- get_number_types_implied_by_object(x, best_types)
    length_x <- g$length
    best_types <- g$types
    if(!is.null(x) & length_x != 0) {
      if(best_guess != 0) {
        # We are looking for the full number of types, it is the largest of the two.
        # It is assumed that if length_x < best_guess, then we will later subset to
        # only some of the types for this object.
        best_guess <- max(best_guess, length_x)
      } else {
        best_guess <- length_x
      }
    }
  }

  if(best_guess == 0) { # Nothing allowed us to guess the number of types, so use default
    list(number_types = default_number_types,
         types = best_types)
  } else {
    list(number_types = best_guess,
         types = best_types)
  }
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
  g <- get_number_types_and_check_conformance(default_number_types = default_number_types,
                                              types = types,
                                              alpha,
                                              gamma,
                                              beta0,
                                              beta,
                                              short_range,
                                              medium_range,
                                              long_range,
                                              ...)

  number_types <- g$number_types
  full_types <- g$types

  if(is.null(full_types)) {
    full_types <- make_default_types(size = number_types)
  }

  alpha <- lapply(alpha, function(a) {
    a <- construct_if_missing(a, 0, number_types, matrix = TRUE)
    if(!isSymmetric(a)) {
      stop("One of the alphas is not symmetric.")
    }
    if(is.null(colnames(a)) & is.null(rownames(a))) {
      colnames(a) <- rownames(a) <- full_types
    }
    a
  })

  gamma <- construct_if_missing(gamma, 0, number_types, matrix = TRUE)
  if(!isSymmetric(gamma)) {
    stop("Gamma is not symmetric.")
  }
  if(is.null(colnames(gamma)) & is.null(rownames(gamma))) {
    colnames(gamma) <- rownames(gamma) <- full_types
  }

  beta0 <- construct_if_missing(beta0, 0, number_types, matrix = FALSE)
  if(is.null(names(beta0))) {
    names(beta0) <- full_types
  }

  beta_nrows <- number_types
  beta_ncols <- length(covariates)

  beta <- construct_if_missing(beta, 1., beta_nrows, beta_ncols, matrix = TRUE)
  if(beta_ncols != 0 && (nrow(beta) != beta_nrows || ncol(beta) != beta_ncols)) {
    print(beta)
    stop("Beta (printed above) does not have the right dimensions.")
  }
  if(is.null(rownames(beta))) {
    rownames(beta) <- full_types
  }

  default_short_distances <- seq(from = 0, to = 0.1, length.out = length(short_range) + 1)[-1]
  short_range <- setNames(lapply(seq_len(length(short_range)), function(i) {
    s <- construct_if_missing(short_range[[i]], default_short_distances[i], number_types, matrix = TRUE)
    if(!isSymmetric(s)) {
      stop("One of the short-range interaction radii matrices is not symmetric.")
    }
    if(is.null(colnames(s)) & is.null(rownames(s))) {
      colnames(s) <- rownames(s) <- full_types
    }
    s
  }), nm = names(short_range))

  medium_range <- construct_if_missing(medium_range, 0, number_types, matrix = TRUE)
  if(!isSymmetric(medium_range)) {
    stop("The medium-range interaction radii matrix is not symmetric.");
  }
  if(is.null(colnames(medium_range)) & is.null(rownames(medium_range))) {
    colnames(medium_range) <- rownames(medium_range) <- full_types
  }

  long_range <- construct_if_missing(long_range, 0, number_types, matrix = TRUE)
  if(!isSymmetric(long_range)) {
    stop("The long-range interaction radii matrix is not symmetric.")
  }
  if(is.null(colnames(long_range)) & is.null(rownames(long_range))) {
    colnames(long_range) <- rownames(long_range) <- full_types
  }

  list(alpha = alpha,
       beta0 = beta0,
       covariates = covariates,
       beta = beta,
       gamma = gamma,
       short_range = short_range,
       medium_range = medium_range,
       long_range = long_range,
       full_types = full_types)
}
#' Generate parameters for the model.
#'
#' Use whatever provided subset of parameters to define defaults for the others.
#' The main unknowns are the intended number of potentials, and the number of types.
#' The intended number of potentials is deduced from the length of alpha, short_range and/or model
#' while the intended number of types is deduced from the dimensions of beta0, alpha, or other parameters.
#' Once these unknowns are found, sensible defaults are used for all parameters..
#'
#' @param window Observation window. Preferably a `ppjsdm` Window, such as `ppjsdm::Rectangle_window`, but also accepts `spatstat` `im` or `owin` objects.
#' @param alpha Matrix of short-range interaction coefficients.
#' Default is a square matrix, filled with zeroes.
#' @param gamma Matrix of medium-range interaction coefficients.
#' Default is a square matrix, filled with zeroes.
#' @param beta0 Vector representing the intercept.
#' Default is a vector filled with zeroes.
#' @param covariates List of covariates. These are converted to the `im` format by applying `as.im` to all elements in the list.
#' Default is an empty list.
#' @param beta Fitted regression coefficients with respect to covariates.
#' Default is a matrix of zeroes, with as many rows as there are types, and number of columns equal to the number of covariates.
#' @param short_range Symmetric matrix of short-range interaction radii. Can also be a list of matrices, each entry representing a different potential.
#' Filled with 0.1 by default, or values equidtsirbuted between 0 and 0.1 if the user intends multiple potentials.
#' @param medium_range Symmetric matrix of medium-range interaction radii. Filled with 0 by default.
#' @param long_range Symmetric matrix of long-range interaction radii. Filled with 0 by default.
#' @param saturation Saturation parameter of the point process. Default is 2.
#' @param types Character vector, with entry i representing the name of type i. Default is a vector (type1, type2, ...) of length the number of types.
#' @param model String representing the short-range model to use. The currently authorised models are obtained with a call to `show_short_range_models()`.
#' @param medium_range_model String representing the medium-range model to use. The currently authorised models are obtained with a call to `show_medium_range_models()`.
#' @param default_number_types Default number of types. If no other over-riding number of types can be deduced from the other parameters, this will be used.
#' @param relevant_types (Optional) Some types that we are interested in. All the model parameters will then be subsetted to include only these types.
#' @param relevant_indices (Optional) Indices of some of the types that we are interested in. All the model parameters will then be subsetted to include only types with these indices.
#' @param ... Other parameters used to infer the number of types. Typically other relevant vectors/matrices that the user has supplied, and could help identify the number of types.
#' @importFrom methods is
#' @importFrom spatstat.geom is.im
#' @examples
#' # Define some of the parameters.
#'
#' # Two short-range potentials:
#' # one with alpha coefficient 0, the second with a 3 x 3 matrix coefficient with 0.1 values.
#' alpha <- list(0, matrix(0.1, 3, 3))
#'
#' # Potentials with an exponential shape.
#' model <- "exponential"
#'
#' # Call the function
#'
#' parameters <- ppjsdm::model_parameters(alpha = alpha, model = model)
#'
#' # Inspect the result
#'
#' # The function has converted alpha into a list of matrices,
#' # deducing that the first entry was also intended to be a matrix.
#' print(parameters$alpha)
#'
#' # The single provided potential shape is assumed to model the shape of both potentials.
#' print(parameters$model)
#'
#' # The rest of the parameters are initialised given the provided information,
#' # so e.g., gamma is initialised as a 3x3 matrix.
#' print(parameters$gamma)
#'
#' @export
#' @md
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
                             relevant_types = NULL,
                             relevant_indices = NULL,
                             ...) {
  # We start with setting all the default arguments.
  # The reason we do not set these as defaults in the function definition
  # is that we want to propagate missing arguments. Indeed, this function will
  # be called from within a range of functions in ppjsdm, and we want to properly detect missing arguments.
  # The other option is to put NULL default arguments everywhere. This seems cleaner to me.
  if(missing(window)) {
    window <- ppjsdm::Rectangle_window()
  }

  if(missing(covariates)) {
    covariates <- list()
  }

  if(missing(saturation)) {
    saturation <- 2
  }

  if(missing(medium_range_model)) {
    medium_range_model <- "square_exponential"
  }

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

  if(missing(model)) {
    model <- NULL
  }

  # At this point, we are done with default parameters.
  # Try to figure out how many short_range potentials the user intends to use
  number_short_range <- 1 # By default, use 1
  if(!is.null(alpha)) {
    # If alpha is not a matrix, then it must be a vector or list of interaction coefficients
    # and we may assume that the length is equal to the intended number of potentials
    if(!is.matrix(alpha)) {
      number_short_range <- length(alpha)
    }
  }
  if(!is.null(short_range)) {
    # Same reasoning as for alpha
    if(!is.matrix(short_range)) {
      if(number_short_range > 1 & length(short_range) > 1 & number_short_range != length(short_range)) {
        stop("Alpha indicated >1 potentials, but short_range does not point to the same number of potentials.")
      }
      number_short_range <- max(number_short_range, length(short_range))
    }
  }
  if(!is.null(model)) {
    # Same reasoning as above
    if(number_short_range > 1 & length(model) > 1 & number_short_range != length(model)) {
      stop("Either alpha or short_range indicated >1 potentials, but model does not point to the same number of potentials.")
    }
    number_short_range <- max(number_short_range, length(model))
  }

  # At this point, we know how many short-range potentials the user intends. We can therefore default-initialise
  # the parameters that depend on that number.
  if(is.null(alpha) | is.matrix(alpha) | (is.numeric(alpha) & length(alpha) == 1)) {
    alpha <- lapply(seq_len(number_short_range), function(x) alpha)
  }

  if(is.null(short_range) | is.matrix(short_range) | (is.numeric(short_range) & length(short_range) == 1)) {
    short_range <- lapply(seq_len(number_short_range), function(x) short_range)
  }

  if(is.null(model) | length(model) == 1 && !is.list(model)) {
    model <- lapply(seq_len(number_short_range), function(x) model)
  }

  # Format some of the types if they are provided
  if(!is.null(alpha)) {
    alpha <- as.list(alpha)
  }
  if(!is.null(short_range)) {
    short_range <- as.list(short_range)
  }
  if(!is.null(model)) {
    model <- as.list(model)
  }
  if(!is.null(beta)) {
    beta <- as.matrix(beta)
  }

  if(!is.null(types)) {
    types <- unlist(types)
  }

  # Model is a list of length the number of potentials, missing values indicated by NULL
  model <- setNames(lapply(model, function(x) {
    if(is.null(x)) {
      "square_bump"
    } else {
      x
    }
  }), nm = names(model))

  # Force window to the right class
  window <- as.Window(window)

  # Make covariates im objects with proper names.
  covariates <- coerce_to_named_im_objects(covariates, "unnamed_covariate", window)

  parameters <- make_default_model_parameters(alpha = alpha,
                                              beta0 = beta0,
                                              covariates = covariates,
                                              beta = beta,
                                              gamma = gamma,
                                              short_range = short_range,
                                              medium_range = medium_range,
                                              long_range = long_range,
                                              types = types,
                                              default_number_types = default_number_types,
                                              ...)

  # Make sure the parameters involving the number of short-range potentials were correctly constructed
  if(length(model) != length(parameters$alpha) || length(model) != length(parameters$short_range)) {
    stop(paste0("Some of the parameters have incompatible sizes; model: [", paste0(model, collapse = ", "),
                "], alpha: [", paste0(parameters$alpha, collapse = ", "),
                "] and short_range: [", paste0(parameters$short_range, collapse = ", "), "]."))
  }

  # See if there exists some potential names
  potential_names <- if(!is.null(names(model))) {
    names(model)
  } else if(!is.null(names(parameters$short_range))) {
    names(parameters$short_range)
  } else {
    NULL
  }
  parameters$potential_names <- potential_names

  if(!is.null(relevant_indices)) {
    relevant <- sort(parameters$full_types[relevant_indices])
  } else {
    relevant <- NULL
  }

  if(!is.null(relevant_types)) {
    relevant <- sort(unique(c(relevant, relevant_types)))
  }

  if(!is.null(relevant)) {
    # If relevant types are supplied, make sure everything is coherent
    # In this case, relevant have to be a subset of the types
    if(all(relevant %in% parameters$full_types)) {
      parameters$alpha <- lapply(parameters$alpha, function(a) {
        z <- as.matrix(a[relevant, relevant])
        colnames(z) <- rownames(z) <- relevant
        z
      })

      parameters$gamma <- as.matrix(parameters$gamma[relevant, relevant])
      colnames(parameters$gamma) <- rownames(parameters$gamma) <- relevant

      parameters$short_range <- lapply(parameters$short_range, function(s) {
        z <- as.matrix(s[relevant, relevant])
        colnames(z) <- rownames(z) <- relevant
        z
      })

      parameters$medium_range <- as.matrix(parameters$medium_range[relevant, relevant])
      colnames(parameters$medium_range) <- rownames(parameters$medium_range) <- relevant

      parameters$long_range <- as.matrix(parameters$long_range[relevant, relevant])
      colnames(parameters$long_range) <- rownames(parameters$long_range) <- relevant

      parameters$beta0 <- parameters$beta0[relevant]

      parameters$beta <- if(ncol(parameters$beta) > 0) {
        if(length(relevant) == 1) { # This avoids some annoying bugs when only one relevant type
          z <- matrix(parameters$beta[relevant, ], nrow = length(relevant))
          colnames(z) <- colnames(parameters$beta)
          z
        } else {
          as.matrix(parameters$beta[relevant, ])
        }
      } else { # And this avoids some annoying bugs when no covariates
        matrix(NA, ncol = 0, nrow = length(relevant))
      }
      rownames(parameters$beta) <- relevant
    } else {
      stop(paste0("The types of the configuration are not a subset of those given by the parameters, configuration: ",
                  paste0(relevant, collapse = ", "),
                  " and supplied types: ",
                  paste0(parameters$full_types, collapse = ", ")))
    }

    # At this point, the parameters should exactly correspond to the types of the configuration
    if(!all(names(parameters$beta0) == relevant)) {
      stop(paste0("Unexpected setting when computing the Papangelou intensity: the supplied parameters do not refer to the same types as the configuration, configuration: ",
                  paste0(relevant, collapse = ", "),
                  " and supplied parameters: ",
                  paste0(names(parameters$beta0), collapse = ", ")))
    }
  }

  # Add the other parameters to the return object
  parameters$window <- window
  parameters$covariates <- covariates
  parameters$saturation <- saturation
  parameters$model <- model
  parameters$medium_range_model <- medium_range_model

  parameters
}
