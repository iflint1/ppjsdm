#' Calculate Variance-Covariance Matrix for a Fitted Gibbs point process.
#'
#' @param object A fitted model object.
#' @param ... Ignored.
#' @export
vcov.gibbsm <- function(object, ...) {
  shortened_regressors <- object$data_list$regressors

  if(!object$estimate_alpha) {
    shortened_regressors <- as.matrix(shortened_regressors[, -grep("alpha_", colnames(shortened_regressors))])
  }
  if(!object$estimate_gamma) {
    shortened_regressors <- as.matrix(shortened_regressors[, -grep("gamma_", colnames(shortened_regressors))])
  }
  if(ncol(shortened_regressors) == 1) { # In this case, R deletes the column names...
    colnames(shortened_regressors) <- "log_lambda1"
  }

  #TODO: Should enforce sizeof(configuration_list) == 1
  #TODO: Should enforce !glmnet
  #TODO: Should enforce constant rho
  rho <- exp(-object$data_list$shift[1])
  fits_coefficients_vector <- object$coefficients_vector
  fits_coefficients <- object$coefficients

  vc <- compute_vcov(configuration = object$configuration_list[[1]],
                       covariates = object$parameters$covariates,
                       model = object$parameters$model,
                       medium_range_model = object$parameters$medium_range_model,
                       short_range = fits_coefficients$short_range,
                       medium_range = fits_coefficients$medium_range,
                       long_range = fits_coefficients$long_range,
                       saturation = object$parameters$saturation,
                       alpha = fits_coefficients$alpha,
                       beta0 = fits_coefficients$beta0,
                       beta = fits_coefficients$beta,
                       gamma = fits_coefficients$gamma,
                       rho = rho,
                       coefficients_vector = fits_coefficients_vector,
                       shortened_regressors = shortened_regressors,
                       data_list = object$data_list,
                       estimate_alpha = object$estimate_alpha,
                       estimate_gamma = object$estimate_gamma)
  vc
}
