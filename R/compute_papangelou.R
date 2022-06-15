#' Compute the Papangelou conditional intensity of the model.
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
#' @param nthreads Maximum number of threads for parallel computing.
#' @export
#' @method compute_papangelou default
compute_papangelou.default <- function(configuration,
                                       x,
                                       y,
                                       type = rep(1, length(x)),
                                       mark = rep(1.0, length(x)),
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
                                       saturation,
                                       nthreads = 1) {

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

  check_type_mark <- function(obj) {
    if(length(obj) != length(x)) {
      if(length(obj) == 1) {
        rep(obj, length(x))
      } else {
        stop("Unknown format for type or mark.")
      }
    } else {
      obj
    }
  }

  type <- check_type_mark(type)
  mark <- check_type_mark(mark)

  compute_papangelou_cpp(x = x,
                         y = y,
                         type = type,
                         mark = mark,
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
                         saturation = parameters$saturation,
                         nthreads = nthreads)
}

#' Compute the Papangelou conditional intensity of the model from a fit object.
#'
#' @param fit Fit object obtained by running gibbsm.
#' @param type Type of the point at which to evaluate the Papangelou conditional intensity.
#' @param configuration Configuration of points at which to evaluate the Papangelou conditional intensity.
#' @param alpha Short range repulsion parameter.
#' @param gamma Medium range repulsion parameter.
#' @param beta0 Log-intensity.
#' @param covariates Covariates.
#' @param beta Covariates coefficients.
#' @param short_range Matrix of short range distances.
#' @param medium_range Matrix of medium range distances.
#' @param long_range Matrix of long range distances.
#' @param saturation Saturation.
#' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
#' @param model Model for short range interaction.
#' @param medium_range_model Model for medium range interaction.
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
