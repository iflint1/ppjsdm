fit_gibbs <- function(gibbsm_data_list, use_glmnet, use_aic) {
  if(use_glmnet) {
    lapply(gibbsm_data_list, function(gibbsm_data) {
      regressors <- gibbsm_data$regressors
      nregressors <- ncol(regressors)
      pfactor <- rep(1, nregressors)
      pfactor[startsWith(colnames(regressors), "shifted_log_lambda")] <- 0

      # Avoid a bug in glmnet: if intercept = FALSE, and there's a column of one  s, it gets ignored by glmnet
      # even though its penalty factor is zero.
      if(all(1 == regressors[, startsWith(colnames(regressors), "shifted_log_lambda")])) {
        regressors[1, startsWith(colnames(regressors), "shifted_log_lambda")] <- 1.001
      }

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

      if(use_aic) {
        coef <- coefficients(fit)[, which.min(aic)]
      } else {
        coef <- coefficients(fit)[, which.min(bic)]
      }

      # We don't use an offset explicitely because the call to glmnet above returns nonsensical results or hangs.
      # Instead, use a shift for all the log_lambda regressors according to -log(rho).
      shift <- gibbsm_data$shift
      number_types <- length(shift)
      lambda <- vector(mode = "numeric", length = number_types)
      alpha <- matrix(0, number_types, number_types)
      gamma <- matrix(0, number_types, number_types)
      for(i in seq_len(number_types)) {
        # Shift all columns with row name shifted_log_lambdai
        lambda[i] <- exp(coef[match(paste0("shifted_log_lambda", i), names(coef))] - shift[i])
        alpha[i, i] <- coef[match(paste0("alpha_", i, "_", i), names(coef))]
        gamma[i, i] <- coef[match(paste0("gamma_", i, "_", i), names(coef))]
        if(i < number_types) {
          for(j in (i + 1):number_types) {
            alpha[i, j] <-  alpha[j, i] <- coef[match(paste0("alpha_", i, "_", j), names(coef))]
            gamma[i, j] <- gamma[j, i] <- coef[match(paste0("gamma_", i, "_", j), names(coef))]
          }
        }
      }

      if(use_aic) {
        list(fit = fit, coefficients = list(lambda = lambda, alpha = alpha, gamma = gamma), aic = min(aic), bic = bic[which.min(aic)])
      } else {
        list(fit = fit, coefficients = list(lambda = lambda, alpha = alpha, gamma = gamma), aic = aic[which.min(bic)], bic = min(bic))
      }
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
