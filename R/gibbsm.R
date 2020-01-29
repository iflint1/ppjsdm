#' Fit a multivariate Gibbs model to a dataset.
#'
#' @param configuration_list A single configuration or a list of configurations assumed to be drawn from the multivariate Gibbs.
#' @param window Observation window.
#' @param covariates Environmental covariates driving the intensity.
#' @param traits Species' traits.
#' @param model String to represent the model we're calibrating. You can check the currently authorised models with a call to `show_model()`.
#' @param radius Interaction radius.
#' @param print Print the fitted coefficients?
#' @param use_glmnet Use `glmnet` instead of `glm`?
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats as.formula binomial coefficients glm.fit lm
#' @importFrom spatstat as.im as.owin
#' @export
gibbsm <- function(configuration_list, window = Rectangle_window(), covariates = list(), traits = list(), model = "identity", radius = NULL, print = TRUE, use_glmnet = TRUE) {
  # Make covariates im objects with proper names.
  covariates <- coerce_to_named_im_objects(covariates, "unnamed_covariate", window)

  # If we're given a single configuration, convert it to a list.
  if(inherits(configuration_list, "Configuration")) {
    configuration_list <- list(configuration_list)
  }

  # Make sure we're given a list of configurations.
  stopifnot(inherits(configuration_list[[1]], "Configuration"))

  gibbsm_data_list <- lapply(configuration_list, function(configuration) {
    prepare_gibbsm_data(configuration, window, covariates, model, radius)
  })
  if(use_glmnet) {
    fits <- lapply(gibbsm_data_list, function(gibbsm_data) {
      glmnet(x = gibbsm_data$regressors,
             y = gibbsm_data$response,
             offset = gibbsm_data$offset,
             alpha = 0.5, # For some reason alpha > 0.9 gives bad results.
             intercept = FALSE,
             family = "binomial")
    })
    fits_coefficients <- lapply(fits, function(fit) coefficients(fit, s = 0))
    cv_fits <- lapply(gibbsm_data_list, function(gibbsm_data) {
      cv.glmnet(x = gibbsm_data$regressors,
                y = gibbsm_data$response,
                offset = gibbsm_data$offset)
    })
  } else {
    fits <- lapply(gibbsm_data_list, function(gibbsm_data) {
      g <- glm.fit(x = gibbsm_data$regressors,
                   y = gibbsm_data$response,
                   offset = gibbsm_data$offset,
                   intercept = FALSE,
                   family = binomial())
      # Note: glm.fit does not correctly set the class,
      # so the user cannot use `glm` methods...
      class(g) <- c(g$class, c("glm", "lm"))
      g
    })
    fits_coefficients <- lapply(fits, function(fit) coefficients(fit))
  }

  number_configurations <- length(configuration_list)

  # If traits have been given, proceed to fit our submodel for alpha ~ traits.
  number_traits <- length(traits)
  if(number_traits > 0) {
    # TODO: Move to C++?
    number_of_individuals_by_configuration <- lapply(configuration_list, function(configuration) get_number_points(configuration, total = FALSE))
    species_names_by_configuration <- lapply(number_of_individuals_by_configuration, function(n) names(n))
    number_of_species_by_configuration <- lapply(number_of_individuals_by_configuration, function(n) length(n))
    number_of_alpha_coefficients_by_configuration <- lapply(number_of_species_by_configuration, function(n) n * (n + 1) / 2)
    regression_length <- Reduce("+", number_of_alpha_coefficients_by_configuration)
    response <- vector(mode = "numeric", length = regression_length)
    joint_regressors <- matrix(data = NA, nrow = regression_length, ncol = number_traits)

    index <- 1
    for(i in seq_len(number_configurations)) {
      for(j in seq_len(number_of_species_by_configuration[[i]])) {
        current_name_j <- species_names_by_configuration[[i]][j]
        for(k in j:number_of_species_by_configuration[[i]]) {
          alpha_string <- paste0("alpha", "_", j, "_", k)
          if(use_glmnet) {
            indices <- which(rownames(fits_coefficients[[i]]) %in% alpha_string)
          } else {
            indices <- which(names(fits_coefficients[[i]]) %in% alpha_string)
          }
          response[index] <- fits_coefficients[[i]][indices]
          current_name_k <- species_names_by_configuration[[i]][k]
          joint_regressors[index, ] <- sapply(traits, function(trait) {
            (trait[current_name_j] - mean(trait)) * (trait[current_name_k] - mean(trait)) / (mean(trait) * mean(trait))
          })
          index <- index + 1
        }
      }
    }

    colnames(joint_regressors) <- paste0("joint_", names(traits))
    formula <- paste0("response ~ 1 + ", paste(colnames(joint_regressors), collapse = " + "))
    traits_fit <- lm(as.formula(formula), data = data.frame(response, joint_regressors))
  }

  if(print) {
    lapply(fits_coefficients, function(fit) print(fit))
    if(number_traits > 0) {
      print(coefficients(traits_fit))
    }
  }

  if(number_configurations == 1) {
    fits <- fits[[1]]
    if(use_glmnet) {
      cv_fits <- cv_fits[[1]]
    }
  }
  ret <- list(complete = fits)
  if(number_traits > 0) {
    ret <- append(list, list(traits = traits_fit))
  }
  if(use_glmnet) {
    ret <- append(ret, list(cv = cv_fits))
  }
  if(length(ret) == 1) {
    ret[[1]]
  } else {
    ret
  }
}
