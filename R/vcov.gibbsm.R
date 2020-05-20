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
  saturation <- object$parameters$saturation

  fits_coefficients <- object$coefficients

  if(!object$estimate_alpha) {
    regressors <- as.matrix(regressors[, -grep("alpha_", colnames(regressors))])
  }
  if(!object$estimate_gamma) {
    regressors <- as.matrix(regressors[, -grep("gamma_", colnames(regressors))])
  }
  if(ncol(regressors) == 1) { # In this case, R deletes the column names...
    colnames(regressors) <- "log_lambda1"
  }

  ntypes <- length(object$data_list$shift)
  #TODO: Should enforce sizeof(configuration_list) == 1
  #TODO: Should enforce constant rho
  #TODO: Is this the correct rho or should we multiply by ntypes?
  rho <- exp(-object$data_list$shift[1])

  fits_coefficients_vector <- object$coefficients_vector
  papangelou <- apply(regressors, 1, function(row) exp(fits_coefficients_vector %*% row))

  t_over_papangelou <- lapply(1:nrow(regressors), function(row) {
    val <- regressors[row, ] / (papangelou[row] + rho)
    beta0 <- val[startsWith(names(val), "log_lambda")]
    if(!object$estimate_alpha) {
      alpha <- rep(0., ntypes * (ntypes + 1) / 2)
    } else {
      alpha <- val[startsWith(names(val), "alpha_")]
    }
    if(!object$estimate_gamma) {
      gamma <- rep(0., ntypes * (ntypes + 1) / 2)
    } else {
      gamma <- val[startsWith(names(val), "gamma_")]
    }
    beta <- val[!(startsWith(names(val), "log_lambda") | startsWith(names(val), "alpha_") | startsWith(names(val), "gamma_"))]
    val <- c(beta0, alpha, gamma, beta)

    list(x = object$data_list$x[row],
         y = object$data_list$y[row],
         type = object$data_list$type[row],
         mark = object$data_list$mark[row],
         value = val)
  })

  vc <- compute_A2_A3(configuration = configuration_list[[1]],
                       covariates = covariates,
                       model = model,
                       medium_range_model = medium_range_model,
                       short_range = fits_coefficients$short_range,
                       medium_range = fits_coefficients$medium_range,
                       long_range = fits_coefficients$long_range,
                       saturation = saturation,
                       alpha = fits_coefficients$alpha,
                       beta0 = fits_coefficients$beta0,
                       beta = fits_coefficients$beta,
                       gamma = fits_coefficients$gamma,
                       rho = rho,
                       t_over_papangelou = t_over_papangelou,
                       coefficients_vector = fits_coefficients_vector,
                       regressors = regressors)

  if(!object$estimate_alpha) {
    keep <- -grep("alpha", colnames(vc))
    vc <- vc[keep, keep]
  }
  if(!object$estimate_gamma) {
    keep <- -grep("gamma", colnames(vc))
    vc <- vc[keep, keep]
  }
  vc
}
