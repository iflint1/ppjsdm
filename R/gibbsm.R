
add_names <- function(str, covariates) {
  if(is.null(names(covariates))) {
    no_name <- rep(TRUE, length(covariates))
  } else {
    no_name <- names(covariates) == ""
  }
  names(covariates)[no_name] <- sprintf(paste0(str, "%d"), seq_len(length(which(no_name))))
  covariates
}

coerce_to_im <- function(lst, window) {
  lapply(as.list(lst), function(element) as.im(element, W = as.owin(window)))
}

#' Fit a multivariate Gibbs model to a dataset.
#'
#' @param configuration_list A single configuration or a list of configurations assumed to be drawn from the multivariate Gibbs.
#' @param window Observation window.
#' @param covariates Environmental covariates driving the intensity.
#' @param traits Species' traits.
#' @param model String to represent the model we're calibrating. You can check the currently authorised models with a call to `show_model()`.
#' @param radius Interaction radius.
#' @param print Print the fitting summary?
#' @importFrom stats as.formula binomial coefficients glm lm
#' @importFrom spatstat as.im as.owin
#' @export
gibbsm <- function(configuration_list, window = Rectangle_window(), covariates = list(), traits = list(), model = "identity", radius = NULL, print = TRUE) {
  covariates <- coerce_to_im(covariates, window)
  covariates <- add_names("covariates", covariates)

  # If we're given a single configuration, convert it to a list.
  if(inherits(configuration_list, "Configuration")) {
    configuration_list <- list(configuration_list)
  }
  # Make sure we're given a list of configurations.
  stopifnot(inherits(configuration_list[[1]], "Configuration"))

  gibbsm_data_list <- lapply(configuration_list, function(configuration) prepare_gibbsm_data(configuration, window, covariates, model, radius))
  glm_fits <- lapply(gibbsm_data_list, function(gibbsm_data) glm(as.formula(gibbsm_data$formula), data = as.data.frame(gibbsm_data$data), family = binomial()))

  number_of_individuals_by_configuration <- lapply(configuration_list, function(configuration) get_number_points(configuration, total = FALSE))
  species_names_by_configuration <- lapply(number_of_individuals_by_configuration, function(n) names(n))
  number_of_species_by_configuration <- lapply(number_of_individuals_by_configuration, function(n) length(n))
  number_configurations <- length(configuration_list)
  number_traits <- length(traits)
  regression_length <- Reduce("+", number_of_species_by_configuration)
  response <- vector(mode = "numeric", length = regression_length)
  regressors <- matrix(data = NA, nrow = regression_length, ncol = number_traits)

  # TODO: Vectorise / C++
  index <- 1
  for(i in seq_len(number_configurations)) {
    for(j in seq_len(number_of_species_by_configuration[[i]])) {
      alpha_string <- paste0("alpha", "_", j, "_", j)
      response[index] <- coefficients(glm_fits[[i]])[alpha_string]
      for(k in seq_len(number_traits)) {
        current_name <- species_names_by_configuration[[i]][j]
        # TODO: Why not access with [current_name] and use vectors?
        val <- traits[[k]][[current_name]]
        regressors[index, k] <- val
      }
      index <- index + 1
    }
  }

  if(print) {
    lapply(glm_fits, function(fit) print(summary(fit)))
    if(number_traits > 0) {
      # TODO: Better names for regressors
      g <- lm(response ~ 1 + regressors)
      print(summary(g))
    }
  }

  if(number_configurations > 1) {
    glm_fits
  } else {
    glm_fits[[1]]
  }
}
