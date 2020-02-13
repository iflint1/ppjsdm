#' Fit a multivariate Gibbs model to a dataset.
#'
#' @param configuration_list A single configuration or a list of configurations assumed to be drawn from the multivariate Gibbs.
#' @param window Observation window.
#' @param covariates Environmental covariates driving the intensity.
#' @param traits Species' traits.
#' @param model String to represent the model we're calibrating. You can check the currently authorised models with a call to `show_short_range_models()`.
#' @param medium_range_model String to represent the model we're calibrating. You can check the currently authorised models with a call to `show_medium_range_models()`.
#' @param short_range Short range interaction radius.
#' @param medium_range Medium range interaction radius.
#' @param long_range Long range interaction radius.
#' @param saturation Saturation parameter of the point process.
#' @param print Print the fitted coefficients?
#' @param use_glmnet Use `glmnet` instead of `glm`?
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats as.formula binomial coefficients glm.fit lm
#' @importFrom spatstat as.im as.owin
#' @export
gibbsm <- function(configuration_list, window = Rectangle_window(), covariates = list(), traits = list(), model = "square_bump", medium_range_model = "square_exponential", short_range = NULL, medium_range = NULL, long_range = NULL, saturation = 2, print = TRUE, use_glmnet = TRUE) {
  # Make covariates im objects with proper names.
  covariates <- coerce_to_named_im_objects(covariates, "unnamed_covariate", window)

  # TODO: This is giving really unexpected results when sizeof(radius) > sizeof(one of the configurations).
  # See the traits part of the fitting vignette. Perhaps add size check for radius?

  # If we're given a single configuration, convert it to a list.
  if(inherits(configuration_list, "Configuration")) {
    configuration_list <- list(configuration_list)
  }

  # Make sure we're given a list of configurations.
  stopifnot(inherits(configuration_list[[1]], "Configuration"))

  gibbsm_data_list <- lapply(configuration_list, function(configuration) {
    prepare_gibbsm_data(configuration, window, covariates, model, medium_range_model, short_range, medium_range, long_range, saturation)
  })
  fitted <- fit_gibbs(gibbsm_data_list, use_glmnet)
  fits <-  lapply(fitted, function(fit) fit$fit)
  fits_coefficients <-lapply(fitted, function(fit) fit$coefficients)
  cv_fits <- lapply(fitted, function(fit) fit$cv)

  number_configurations <- length(configuration_list)

  # If traits have been given, proceed to fit our submodel for alpha ~ traits.
  number_traits <- length(traits)
  if(number_traits > 0) {
    # TODO: Move to C++?
    number_of_individuals_by_configuration <- lapply(configuration_list, function(configuration) get_number_points(configuration, total = FALSE))
    species_names_by_configuration <- lapply(number_of_individuals_by_configuration, function(n) names(n))
    number_of_species_by_configuration <- lapply(number_of_individuals_by_configuration, function(n) length(n))

    number_of_nondiagonal_alpha_coefficients_by_configuration <- lapply(number_of_species_by_configuration, function(n) n * (n - 1) / 2)
    number_of_diagonal_alpha_coefficients_by_configuration <- lapply(number_of_species_by_configuration, function(n) n)

    diagonal_regression_length <- Reduce("+", number_of_diagonal_alpha_coefficients_by_configuration)
    nondiagonal_regression_length <- Reduce("+", number_of_nondiagonal_alpha_coefficients_by_configuration)

    diagonal_response <- vector(mode = "numeric", length = diagonal_regression_length)
    nondiagonal_response <- vector(mode = "numeric", length = nondiagonal_regression_length)

    regressors <- matrix(data = NA, nrow = diagonal_regression_length, ncol = number_traits)
    joint_regressors <- matrix(data = NA, nrow = nondiagonal_regression_length, ncol = number_traits)

    normalisation <- sapply(traits, function(trait) max(abs(trait)))
    joint_normalisation <- sapply(traits, function(trait) max(abs(outer(trait, trait, '-'))))

    joint_index <- 1
    index <- 1
    for(i in seq_len(number_configurations)) {
      for(j in seq_len(number_of_species_by_configuration[[i]])) {
        current_name_j <- species_names_by_configuration[[i]][j]

        alpha_string <- paste0("alpha", "_", j, "_", j)
        if(use_glmnet) {
          indices <- which(rownames(fits_coefficients[[i]]) %in% alpha_string)
        } else {
          indices <- which(names(fits_coefficients[[i]]) %in% alpha_string)
        }
        diagonal_response[index] <- fits_coefficients[[i]][indices]
        for(l in seq_len(number_traits)) {
          regressors[index, l] <- traits[[l]][current_name_j] / normalisation[l]
        }
        index <- index + 1

        if(j < number_of_species_by_configuration[[i]]) {
          for(k in (j + 1):number_of_species_by_configuration[[i]]) {
            alpha_string <- paste0("alpha", "_", j, "_", k)
            if(use_glmnet) {
              indices <- which(rownames(fits_coefficients[[i]]) %in% alpha_string)
            } else {
              indices <- which(names(fits_coefficients[[i]]) %in% alpha_string)
            }
            nondiagonal_response[joint_index] <- fits_coefficients[[i]][indices]
            current_name_k <- species_names_by_configuration[[i]][k]
            # TODO: Very inefficient, hopefully move to C++...
            for(l in seq_len(number_traits)) {
              joint_regressors[joint_index, l] <- abs(traits[[l]][current_name_j] - traits[[l]][current_name_k]) / joint_normalisation[l]
            }
            joint_index <- joint_index + 1
          }
        }
      }
    }

    colnames(joint_regressors) <- paste0("joint_", names(traits))
    formula_joint <- paste0("nondiagonal_response ~ 1 + ", paste(colnames(joint_regressors), collapse = " + "))
    joint_traits_fit <- lm(as.formula(formula_joint), data = data.frame(nondiagonal_response, joint_regressors))

    colnames(regressors) <- names(traits)
    formula <- paste0("diagonal_response ~ 1 + ", paste(colnames(regressors), collapse = " + "))
    traits_fit <- lm(as.formula(formula), data = data.frame(diagonal_response, regressors))
  }

  if(print) {
    lapply(fits_coefficients, function(fit) print(fit))
    if(number_traits > 0) {
      print(coefficients(traits_fit))
      print(coefficients(joint_traits_fit))
    }
  }

  if(number_configurations == 1) {
    fits <- fits[[1]]
    fits_coefficients <- fits_coefficients[[1]]
    if(use_glmnet) {
      cv_fits <- cv_fits[[1]]
    }
  }
  ret <- list(complete = fits, coefficients = fits_coefficients)
  if(number_traits > 0) {
    ret <- append(list, list(traits = traits_fit, joint_traits = joint_traits_fit))
  }
  if(use_glmnet) {
    ret <- append(ret, list(cv = cv_fits))
  }
  ret
}
