
#' Compute Papangelou conditional intensity of the model.
#'
#' @param configuration Configuration.
#' @param x Coordinates along the x-axis of the points at which to evaluate the Papangelou conditional intensity.
#' @param y Coordinates along the y-axis of the points at which to evaluate the Papangelou conditional intensity.
#' @param type Type of the point (as an integer >= 1).
#' @param model String representing the model to use. You can check the currently authorised models with a call to `show_models()`.
#' @param medium_range_model String representing the medium range model to use. You can check the currently authorised models with a call to `show_medium_range_models()`.
#' @param alpha Short range repulsion matrix.
#' @param beta0 A vector representing the log-intensities of the point processes.
#' @param beta Coefficients related to covariates.
#' @param gamma Medium range repulsion matrix.
#' @param covariates Covariates.
#' @param short_range Short range interaction radii.
#' @param medium_range Medium range interaction radii.
#' @param long_range Long range interaction radii.
#' @param saturation Saturation parameter.
#' @importFrom stats na.omit
#' @param mark Mark of the point to add.
#' @export
compute_papangelou <- function(configuration,
                               x,
                               y,
                               type,
                               mark,
                               model,
                               medium_range_model,
                               alpha,
                               beta0,
                               beta,
                               gamma,
                               covariates,
                               short_range,
                               medium_range,
                               long_range,
                               saturation) {
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
                                 medium_range_model = medium_range_model)
  # Arguments we want to forward
  fwd_args <- c("type", "mark")
  # Obtain the list of arguments provided
  other_args <- as.list(match.call())
  # Remove first list element, it's the function call
  other_args[[1]] <- NULL
  # Remove the arguments that are not listed in fwd_args
  other_args <- other_args[na.omit(match(fwd_args, names(other_args)))]

  args <- c(other_args, list(x = x,
                             y = y,
                             configuration = as.Configuration(configuration),
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
                             saturation = parameters$saturation))
  do.call(compute_papangelou_cpp, args, envir = parent.frame())
}
