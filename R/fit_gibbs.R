fit_gibbs <- function(gibbsm_data_list, use_glmnet = TRUE) {
  if(use_glmnet) {
    lapply(gibbsm_data_list, function(gibbsm_data) {
      regressors <- gibbsm_data$regressors
      nregressors <- ncol(regressors)
      pfactor <- rep(1, nregressors)
      pfactor[startsWith(colnames(regressors), "shifted_log_lambda")] <- 0

      fit <- glmnet(x = regressors,
                    y = gibbsm_data$response,
                    alpha = 1,
                    intercept = FALSE,
                    family = "binomial",
                    penalty.factor = pfactor)

      tLL <- fit$nulldev - deviance(fit)
      k <- fit$df
      n <- fit$nobs
      aic <- -tLL + 2 * k + 2 * k * (k + 1) / (n - k - 1)
      bic <- log(n) * k - tLL

      coef <- coefficients(fit)[, min(aic) == aic]

      # We don't use an offset explicitely because the call to glmnet above returns nonsensical results or hangs.
      # Instead, use a shift for all the log_lambda regressors according to -log(rho).
      shift <- gibbsm_data$shift
      number_types <- length(shift)
      for(i in seq_len(number_types)) {
        # Shift all columns with row name shifted_log_lambdai
        v <- coef[paste0("shifted_log_lambda", i)] - shift[i]
        coef[paste0("shifted_log_lambda", i)] <- v

        # Remove "shifted" in front of log_lambda names.
        names(coef)[match(paste0("shifted_log_lambda", i), names(coef))] <- paste0("log_lambda", i)
      }

      list(fit = fit, coefficients = coef, aic = min(aic), bic = bic[min(aic) == aic])
    })
  } else {
    lapply(gibbsm_data_list, function(gibbsm_data) {
      fit <- glm.fit(x = gibbsm_data$regressors,
                     y = gibbsm_data$response,
                     offset = gibbsm_data$offset,
                     intercept = FALSE,
                     family = binomial())
      # Note: glm.fit does not correctly set the class,
      # so the user cannot use `glm` methods...
      class(fit) <- c(fit$class, c("glm", "lm"))

      coef <- coefficients(fit)

      list(fit = fit, coefficients = coef)
    })
  }
}
