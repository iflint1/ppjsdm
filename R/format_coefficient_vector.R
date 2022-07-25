#' Format the coefficient vector returned by gibbsm into vector and matrix form for all coefficients.
#'
#' @param coefficient_vector Coefficient vector to format.
#' @param number_types Number of types.
#' @param types_names Names of the types.
#' @param covariates_names Names of the covariates.
#' @param estimate_alpha Matrix where entry i, j is a boolean representing whether or not we are modelling short-range interactions between i and j.
#' @param estimate_gamma Matrix where entry i, j is a boolean representing whether or not we are modelling medium-range interactions between i and j.
#' @export
format_coefficient_vector <- function(coefficient_vector,
                                      number_types,
                                      types_names,
                                      covariates_names,
                                      estimate_alpha,
                                      estimate_gamma) {

  # Construct coefficients in matrix form
  beta0 <- vector(mode = "numeric", length = number_types)
  alpha <- lapply(seq_len(length(estimate_alpha)), function(n) matrix(0, number_types, number_types))
  gamma <- matrix(0, number_types, number_types)

  # Environmental covariates in vector form
  beta_vector <- coefficient_vector[!(startsWith(names(coefficient_vector), "beta0") |
                                        startsWith(names(coefficient_vector), "alpha") |
                                        startsWith(names(coefficient_vector), "gamma") |
                                        ("(Intercept)" == names(coefficient_vector)))]

  # Environmental covariates matrix
  beta <- matrix(NA, nrow = number_types, ncol = length(beta_vector) / number_types)

  # Fill all the coefficients
  for(i in seq_len(number_types)) {
    beta0[i] <- coefficient_vector[match(paste0("beta0_", i), names(coefficient_vector))]
    if(length(beta_vector) != 0) {
      beta[i, ] <- beta_vector[endsWith(names(beta_vector), paste0("_", i))]
    }
    for(k in seq_len(length(estimate_alpha))) {
      if(estimate_alpha[[k]][i, i]) {
        alpha[[k]][i, i] <- coefficient_vector[match(paste0("alpha", k, "_", i, "_", i), names(coefficient_vector))]
      }
    }
    if(estimate_gamma[i, i]) {
      gamma[i, i] <- coefficient_vector[match(paste0("gamma_", i, "_", i), names(coefficient_vector))]
    }
    if(i < number_types) {
      for(j in (i + 1):number_types) {
        for(k in seq_len(length(estimate_alpha))) {
          if(estimate_alpha[[k]][i, j]) {
            alpha[[k]][i, j] <- alpha[[k]][j, i] <- coefficient_vector[match(paste0("alpha", k, "_", i, "_", j), names(coefficient_vector))]
          }
        }
        if(estimate_gamma[i, j]) {
          gamma[i, j] <- gamma[j, i] <- coefficient_vector[match(paste0("gamma_", i, "_", j), names(coefficient_vector))]
        }
      }
    }
  }

  # Fill names
  if(!missing(types_names)) {
    names(beta0) <- types_names
    rownames(beta) <- types_names
    alpha <- lapply(alpha, function(a) {
      z <- a
      colnames(z) <- types_names
      rownames(z) <- types_names
      z
    })
    rownames(gamma) <- types_names
    colnames(gamma) <- types_names
  }
  if(!missing(covariates_names)) {
    colnames(beta) <- covariates_names
  }

  list(beta0 = beta0,
       alpha = alpha,
       gamma = gamma,
       beta = beta)
}
