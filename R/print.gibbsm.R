#' @method print gibbsm
#' @export
print.gibbsm <- function(x, ...) {
  str <- paste0("A fitted saturation pairwise interaction Gibbs point process model.\n\nNumber of replicated samples used for the fit: ",
                length(x$configuration_list),
                ".\nNumber of types: ",
                length(x$coefficients$beta0),
                ".\nNumber of covariates: ",
                length(x$parameters$covariates),
                ".\nShort-range potential(s): ",
                paste0(x$parameters$model, collapse = ", "),
                ".\n")

  if(sum(x$estimate_gamma) > 0) {
    str <- paste0(str, "Using a medium-range potential: ",
                  x$parameters$medium_range_model,
                  ".\n")
  } else {
    str <- paste0(str, "Using no medium-range potential.\n")
  }

  str <- paste0(str, "Saturation parameter: ",
                x$parameters$saturation,
                ".\nWindow: ")

  cat(str)
  print(x$window)
  str <- paste0("Fitting algorithm used: ",
                x$fit_algorithm,
                ".\nAIC of the underlying logistic regression: ",
                x$aic,
                ".\nDistribution of dummy points: ",
                x$dummy_distribution,
                ".\nNumber of dummy points by type: ",
                paste0(paste0(levels(x$data_list$dummy$types), " = ", table(x$data_list$dummy$types)), collapse = ", "),
                ".\n")

  if(x$used_regularization) {
    str <- paste0(str, "With regularisation.\n")
  } else {
    str <- paste0(str, "No regularisation.\n")
  }

  str <- paste0(str, "\n---------------------\nFitted coefficients\n---------------------\n\n",
                "Intercept beta0: \n",
                paste0(paste0(names(x$coefficients$beta0), " = ", as.table(x$coefficients$beta0)), collapse = ", "),
                ".\n\n")

  for(cov in colnames(x$coefficients$beta)) {
    str <- paste0(str, "Regression coefficient with respect to ", cov, ": \n",
                  paste0(paste0(names(x$coefficients$beta0), " = ", as.table(x$coefficients$beta[, cov])), collapse = ", "),
                  ".\n\n")
  }

  cat(str)

  for(i in seq_along(x$coefficients$alpha)) {
    str <- paste0("Short-range interaction coefficients for potential ", i, ": \n")
    cat(str)
    print(x$coefficients$alpha[[i]])
    cat("\n")
  }

  if(sum(x$estimate_gamma) > 0) {
    str <- paste0("Medium-range interaction coefficients : \n")
    cat(str)
    print(x$coefficients$gamma)
  }

}
