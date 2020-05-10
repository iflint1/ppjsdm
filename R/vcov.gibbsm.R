#' Calculate Variance-Covariance Matrix for a Fitted Gibbs point process.
#'
#' @param object A fitted model object.
#' @param ... Ignored.
#' @export
vcov.gibbsm <- function(object, ...) {
  configuration_list <- object$configuration_list

  regressors <- object$data_list$regressors
  model <- object$parameters$model
  medium_range_model <- object$parameters$medium_range_model
  covariates <- object$parameters$covariates
  short_range <- object$parameters$short_range
  medium_range <- object$parameters$medium_range
  long_range <- object$parameters$long_range
  saturation <- object$parameters$saturation

  fits_coefficients <- object$coefficients

  if(!object$estimate_alpha) {
    regressors <- as.matrix(regressors[, -grep("alpha", colnames(regressors))])
  }
  if(!object$estimate_gamma) {
    regressors <- as.matrix(regressors[, -grep("gamma", colnames(regressors))])
  }
  if(ncol(regressors) == 1) { # In this case, R deletes the column names...
    colnames(regressors) <- "shifted_log_lambda1"
  }

  #TODO: Should be different when estimating radii
  #TODO: Should enforce sizeof(configuration_list) == 1
  #TODO: Should enforce constant rho
  rho <- exp(-object$data_list$shift[1])
  S <- rho / window_volume(object$window) * Reduce('+', lapply(1:nrow(regressors), function(row) {
    type_string <- levels(types(configuration_list[[1]]))[object$data_list$type[row]]
    conf <- remove_from_configuration(configuration_list[[1]], c(object$data_list$x[row], object$data_list$y[row]), type_string)
    papangelou <- compute_papangelou(configuration = conf,
                                     x = object$data_list$x[row],
                                     y = object$data_list$y[row],
                                     type = object$data_list$type[row],
                                     model = model,
                                     medium_range_model = medium_range_model,
                                     alpha = fits_coefficients$alpha,
                                     beta0 = fits_coefficients$beta0,
                                     beta = fits_coefficients$beta,
                                     gamma = fits_coefficients$gamma,
                                     covariates = covariates,
                                     short_range = short_range,
                                     medium_range = medium_range,
                                     long_range = long_range,
                                     saturation = saturation,
                                     mark = object$data_list$mark[row])
    regressors[row, ] %*% t(regressors[row, ]) * papangelou / (papangelou + rho)^2
  }))

  A1 <- rho^2 / window_volume(object$window) * Reduce('+', lapply(1:nrow(regressors), function(row) {
    type_string <- levels(types(configuration_list[[1]]))[object$data_list$type[row]]
    conf <- remove_from_configuration(configuration_list[[1]], c(object$data_list$x[row], object$data_list$y[row]), type_string)
    papangelou <- compute_papangelou(configuration = conf,
                                     x = object$data_list$x[row],
                                     y = object$data_list$y[row],
                                     type = object$data_list$type[row],
                                     model = model,
                                     medium_range_model = medium_range_model,
                                     alpha = fits_coefficients$alpha,
                                     beta0 = fits_coefficients$beta0,
                                     beta = fits_coefficients$beta,
                                     gamma = fits_coefficients$gamma,
                                     covariates = covariates,
                                     short_range = short_range,
                                     medium_range = medium_range,
                                     long_range = long_range,
                                     saturation = saturation,
                                     mark = object$data_list$mark[row])
   regressors[row, ] %*% t(regressors[row, ]) * papangelou / (papangelou + rho)^3
  }))

  kappa <- 1 / window_volume(object$window) * Reduce('+', lapply(1:nrow(regressors), function(row) {
    type_string <- levels(types(configuration_list[[1]]))[object$data_list$type[row]]
    conf <- remove_from_configuration(configuration_list[[1]], c(object$data_list$x[row], object$data_list$y[row]), type_string)
    papangelou <- compute_papangelou(configuration = conf,
                                     x = object$data_list$x[row],
                                     y = object$data_list$y[row],
                                     type = object$data_list$type[row],
                                     model = model,
                                     medium_range_model = medium_range_model,
                                     alpha = fits_coefficients$alpha,
                                     beta0 = fits_coefficients$beta0,
                                     beta = fits_coefficients$beta,
                                     gamma = fits_coefficients$gamma,
                                     covariates = covariates,
                                     short_range = short_range,
                                     medium_range = medium_range,
                                     long_range = long_range,
                                     saturation = saturation,
                                     mark = object$data_list$mark[row])
    1. / (papangelou + rho)
  }))

  temp_A1 <- rho^2 / window_volume(object$window) * Reduce('+', lapply(1:nrow(regressors), function(row) {
    type_string <- levels(types(configuration_list[[1]]))[object$data_list$type[row]]
    conf <- remove_from_configuration(configuration_list[[1]], c(object$data_list$x[row], object$data_list$y[row]), type_string)
    papangelou <- compute_papangelou(configuration = conf,
                                     x = object$data_list$x[row],
                                     y = object$data_list$y[row],
                                     type = object$data_list$type[row],
                                     model = model,
                                     medium_range_model = medium_range_model,
                                     alpha = fits_coefficients$alpha,
                                     beta0 = fits_coefficients$beta0,
                                     beta = fits_coefficients$beta,
                                     gamma = fits_coefficients$gamma,
                                     covariates = covariates,
                                     short_range = short_range,
                                     medium_range = medium_range,
                                     long_range = long_range,
                                     saturation = saturation,
                                     mark = object$data_list$mark[row])
    regressors[row, ] %*% t(regressors[row, ]) * papangelou^2 / (papangelou + rho)^3
  }))

  other_temp_A1 <- rho / window_volume(object$window) * Reduce('+', lapply(1:nrow(regressors), function(row) {
    type_string <- levels(types(configuration_list[[1]]))[object$data_list$type[row]]
    conf <- remove_from_configuration(configuration_list[[1]], c(object$data_list$x[row], object$data_list$y[row]), type_string)
    papangelou <- compute_papangelou(configuration = conf,
                                     x = object$data_list$x[row],
                                     y = object$data_list$y[row],
                                     type = object$data_list$type[row],
                                     model = model,
                                     medium_range_model = medium_range_model,
                                     alpha = fits_coefficients$alpha,
                                     beta0 = fits_coefficients$beta0,
                                     beta = fits_coefficients$beta,
                                     gamma = fits_coefficients$gamma,
                                     covariates = covariates,
                                     short_range = short_range,
                                     medium_range = medium_range,
                                     long_range = long_range,
                                     saturation = saturation,
                                     mark = object$data_list$mark[row])
    regressors[row, ] %*% t(rep.int(1, ncol(regressors))) * papangelou / (papangelou + rho)^2
  }))

  t_over_papangelou <- lapply(1:nrow(regressors), function(row) {
    type_string <- levels(types(configuration_list[[1]]))[object$data_list$type[row]]
    conf <- remove_from_configuration(configuration_list[[1]], c(object$data_list$x[row], object$data_list$y[row]), type_string)
    papangelou <- compute_papangelou(configuration = conf,
                                     x = object$data_list$x[row],
                                     y = object$data_list$y[row],
                                     type = object$data_list$type[row],
                                     model = model,
                                     medium_range_model = medium_range_model,
                                     alpha = fits_coefficients$alpha,
                                     beta0 = fits_coefficients$beta0,
                                     beta = fits_coefficients$beta,
                                     gamma = fits_coefficients$gamma,
                                     covariates = covariates,
                                     short_range = short_range,
                                     medium_range = medium_range,
                                     long_range = long_range,
                                     saturation = saturation,
                                     mark = object$data_list$mark[row])
    list(x = object$data_list$x[row],
         y = object$data_list$y[row],
         type = object$data_list$type[row],
         mark = object$data_list$mark[row],
         value = regressors[row, ] / (papangelou + rho))
  })

  A2_A3 <- compute_A2_A3(configuration = configuration_list[[1]],
                         covariates = covariates,
                         model = model,
                         medium_range_model = medium_range_model,
                         short_range = short_range,
                         medium_range = medium_range,
                         long_range = long_range,
                         saturation = saturation,
                         alpha = fits_coefficients$alpha,
                         beta0 = fits_coefficients$beta0,
                         beta = fits_coefficients$beta,
                         gamma = fits_coefficients$gamma,
                         rho = rho,
                         area = window_volume(object$window),
                         t_over_papangelou = t_over_papangelou)
  A2 <- A2_A3$A2
  A3 <- A2_A3$A3
  if(!object$estimate_alpha) {
    keep <- -grep("alpha", colnames(A2))
    A2 <- A2[keep, keep]
    A3 <- A3[keep, keep]
  }
  if(!object$estimate_gamma) {
    keep <- -grep("gamma", colnames(A2))
    A2 <- A2[keep, keep]
    A3 <- A3[keep, keep]
  }

  G2 <- (temp_A1 - other_temp_A1 * t(other_temp_A1) / kappa) / rho
  #G2 <- temp_A1 / rho

  vc <- solve(S) * (A1 + A2 + A3 + G2) * solve(S)
  vc
}
