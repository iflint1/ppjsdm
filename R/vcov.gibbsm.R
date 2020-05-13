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
  #TODO: Should be different when estimating radii
  #TODO: Should enforce sizeof(configuration_list) == 1
  #TODO: Should enforce constant rho
  #TODO: Is this the correct rho or should we multiply by ntypes?
  rho <- exp(-object$data_list$shift[1])

  papangelou <- vector(mode = "numeric", length = nrow(regressors))
  for(row in seq_len(length(papangelou))) {
    type_string <- levels(types(configuration_list[[1]]))[object$data_list$type[row]]
    conf <- remove_from_configuration(configuration_list[[1]], c(object$data_list$x[row], object$data_list$y[row]), type_string)
    papangelou[row] <- compute_papangelou(configuration = conf,
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
                                          short_range = fits_coefficients$short_range,
                                          medium_range = fits_coefficients$medium_range,
                                          long_range = fits_coefficients$long_range,
                                          saturation = saturation,
                                          mark = object$data_list$mark[row])
  }

  S <- Reduce('+', lapply(1:nrow(regressors), function(row) {
    regressors[row, ] %*% t(regressors[row, ]) * papangelou[row] / (papangelou[row] + rho)^2
  }))
  A1 <- Reduce('+', lapply(1:nrow(regressors), function(row) {
    regressors[row, ] %*% t(regressors[row, ]) * papangelou[row] / (papangelou[row] + rho)^3
  }))
  kappa <- sum(1. / (papangelou + rho))
  temp_A1 <- Reduce('+', lapply(1:nrow(regressors), function(row) {
    regressors[row, ] %*% t(regressors[row, ]) * papangelou[row]^2 / (papangelou[row] + rho)^3
  }))

  other_temp_A1 <- Reduce('+', lapply(1:nrow(regressors), function(row) {
    regressors[row, ] %*% t(rep.int(1, ncol(regressors))) * papangelou[row] / (papangelou[row] + rho)^2
  }))

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

  A2_A3 <- compute_A2_A3(configuration = configuration_list[[1]],
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
  # G2 <- temp_A1 / rho

  Sinv <- solve(S)
  vc <- Sinv %*% (A1 + A2 + A3 + G2) %*% Sinv
  # print(S)
  # print(A1)
  # print(A2)
  # print(A3)
  # print(G2)
  vc
}
