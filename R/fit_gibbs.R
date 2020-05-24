fit_gibbs <- function(gibbsm_data, use_glmnet, use_aic, estimate_alpha, estimate_gamma) {
  regressors <- gibbsm_data$regressors

  if(use_glmnet) {
    nregressors <- ncol(regressors)
    pfactor <- rep(1, nregressors)
    pfactor[startsWith(colnames(regressors), "log_lambda")] <- 0

    # Avoid a bug in glmnet: if intercept = FALSE, and there's a column of ones, it gets ignored by glmnet
    # even though its penalty factor is zero.
    if(all(1 == regressors[, startsWith(colnames(regressors), "log_lambda")])) {
      regressors[1, startsWith(colnames(regressors), "log_lambda")] <- 1.001
    }

    # We don't use an offset explicitely because the call to glmnet above returns nonsensical results or hangs.
    # Instead, use a shift for all the log_lambda regressors according to -log(rho).
    fit <- glmnet(x = regressors,
                  y = gibbsm_data$response,
                  alpha = 1,
                  intercept = FALSE,
                  family = "binomial",
                  penalty.factor = pfactor)

    shift <- gibbsm_data$shift

    tLL <- fit$nulldev - deviance(fit)
    k <- fit$df
    n <- fit$nobs
    aic <- -tLL + 2 * k + 2 * k * (k + 1) / (n - k - 1)
    bic <- log(n) * k - tLL

    if(use_aic) {
      coef <- coefficients(fit)[, which.min(aic)]
      aic <- min(aic)
      bic <- bic[which.min(aic)]
    } else {
      coef <- coefficients(fit)[, which.min(bic)]
      aic <- aic[which.min(bic)]
      bic <- min(bic)
    }
    coef <- coef[-which(names(coef) == "(Intercept)")]
    fit_algorithm <- "glmnet"
  } else {
    fit <- glm.fit(x = regressors,
                   y = gibbsm_data$response,
                   offset = gibbsm_data$offset,
                   intercept = FALSE,
                   family = binomial())
    shift <- rep(0, length(gibbsm_data$shift))
    # Note: glm.fit does not correctly set the class,
    # so the user cannot use `glm` methods...
    class(fit) <- c(fit$class, c("glm", "lm"))

    aic <- AIC(fit)
    bic <- BIC(fit)

    coef <- coefficients(fit)
    fit_algorithm <- "glm"
  }

  number_types <- length(shift)
  beta0 <- vector(mode = "numeric", length = number_types)
  alpha <- matrix(0, number_types, number_types)
  gamma <- matrix(0, number_types, number_types)

  beta_vector <- coef[!(startsWith(names(coef), "log_lambda") | startsWith(names(coef), "alpha") | startsWith(names(coef), "gamma") | ("(Intercept)" == names(coef)))]
  if(length(beta_vector) == 0) {
    beta_vector <- matrix(NA, nrow = number_types, ncol = 0)
  }
  beta <- matrix(NA, nrow = number_types, ncol = length(beta_vector) / number_types)
  for(i in seq_len(number_types)) {
    # Shift all columns with row name log_lambdai
    beta0[i] <- coef[match(paste0("log_lambda", i), names(coef))] - shift[i]
    if(length(beta_vector) != 0) {
      beta[i, ] <- beta_vector[endsWith(names(beta_vector), paste0("_", i))]
    }
    if(estimate_alpha) {
      alpha[i, i] <- coef[match(paste0("alpha_", i, "_", i), names(coef))]
    }
    if(estimate_gamma) {
      gamma[i, i] <- coef[match(paste0("gamma_", i, "_", i), names(coef))]
    }
    if(i < number_types) {
      for(j in (i + 1):number_types) {
        if(estimate_alpha) {
          alpha[i, j] <- alpha[j, i] <- coef[match(paste0("alpha_", i, "_", j), names(coef))]
        }
        if(estimate_gamma) {
          gamma[i, j] <- gamma[j, i] <- coef[match(paste0("gamma_", i, "_", j), names(coef))]
        }
      }
    }
  }
  list(fit = fit,
       coefficients = list(beta0 = beta0,
                           alpha = alpha,
                           gamma = gamma,
                           beta = beta),
       coefficients_vector = coef,
       aic = aic,
       bic = bic,
       fit_algorithm = fit_algorithm)
}
