#' Fit a multivariate Gibbs model to a dataset.
#'
#' @param configuration A configuration assumed to be a draw from the multivariate Gibbs.
#' @param window Observation window.
#' @param covariates Environmental covariates driving the intensity.
#' @param traits Species' traits.
#' @param model String to represent the model we're calibrating. You can check the currently authorised models with a call to `show_model()`.
#' @param radius Interaction radius.
#' @param print Print the fitting summary?
#' @importFrom stats as.formula binomial glm
#' @export
gibbsm <- function(configuration, window = Rectangle_window(), covariates = list(), traits = list(), model = "identity", radius = NULL, print = TRUE) {
  ret <- prepare_gibbsm_data(configuration, window, covariates, traits, model, radius)
  g <- glm(as.formula(ret$formula), data = as.data.frame(ret$data), family = binomial())

  if(print) {
    print(summary(g))
  }

  g
}
