fit_gibbs <- function(gibbsm_data_list, use_glmnet = TRUE) {
  if(use_glmnet) {
    fits <- lapply(gibbsm_data_list, function(gibbsm_data) {
      regressors <- gibbsm_data$regressors
      nregressors <- length(regressors)
      pfactor <- rep(1, nregressors)
      pfactor[startsWith(colnames(regressors), "log_lambda")] <- 0

      fit <- glmnet(x = gibbsm_data$regressors,
                    y = gibbsm_data$response,
                    alpha = 1,
                    intercept = FALSE,
                    family = "binomial",
                    penalty.factor = pfactor)
      # We don't use an offset explicitely because the call to glmnet above returns nonsensical results or hangs.
      # Instead, use a shift for all the log_lambda regressors according to -log(rho).
      shift <- gibbsm_data$shift
      number_types <- length(shift)
      for(i in seq_len(number_types)) {
        v <- fit$beta[paste0("log_lambda_", i), ]
        v[v != 0] <- v[v != 0] - shift[i]
        fit$beta[paste0("log_lambda_", i), ] <- v
      }
      fit
    })
    coefficients <- lapply(fits, function(fit) coefficients(fit, s = 0))
    cv <- lapply(gibbsm_data_list, function(gibbsm_data) {
      cv.glmnet(x = gibbsm_data$regressors,
                y = gibbsm_data$response,
                offset = gibbsm_data$offset)
    })
  } else {
    fits <- lapply(gibbsm_data_list, function(gibbsm_data) {
      g <- glm.fit(x = gibbsm_data$regressors,
                   y = gibbsm_data$response,
                   intercept = FALSE,
                   family = binomial())
      # Note: glm.fit does not correctly set the class,
      # so the user cannot use `glm` methods...
      class(g) <- c(g$class, c("glm", "lm"))
      g
    })
    coefficients <- lapply(fits, function(fit) coefficients(fit))
    cv <- list()
  }
  list(fits = fits, coefficients = coefficients, cv = cv)
}
