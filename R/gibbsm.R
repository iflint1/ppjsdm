
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
#' @importFrom stats as.formula binomial glm
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

  if(print) {
    print(summary(glm_fits[[1]]))
  }

  glm_fits[[1]]
}
