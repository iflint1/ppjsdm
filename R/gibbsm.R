#' Fit a multivariate Gibbs model to a dataset.
#'
#' @param configuration A configuration assumed to be a draw from the multivariate Gibbs.
#' @param window Observation window.
#' @param model String to represent the model we're calibrating. At the moment, either "identity", "Strauss" or "Geyer".
#' @param radius Interaction radius.
#' @param print Print the fitting summary?
#' @importFrom stats as.formula binomial glm
#' @export
gibbsm <- function(configuration, window, model = "identity", radius = 0, print = TRUE) {
  # This is the guideline from the Baddeley et al. paper
  rho <- 4 * get_number_points(configuration) / window_volume(window)

  types <- types(configuration)
  distinct_types <- levels(types)
  p <- length(distinct_types)

  # TODO: binomial process should give better results, easy change.
  D <- rppp(window, lambda = rho, types = distinct_types)

  n_Z <- sum(get_number_points(configuration))
  n_D <- sum(get_number_points(D))
  response <- c(rep.int(1, n_Z), rep.int(0, n_D))

  log_lambda <- matrix(0, n_Z + n_D, p)
  # TODO: Horrible, vectorise
  for(i in seq_len(n_Z + n_D)) {
    for(j in seq_len(p)) {
      if(i <= n_Z) {
        if(types(configuration)[i] == distinct_types[j]) {
          log_lambda[i, j] <- 1
        }
      } else {
        if(types(D)[i - n_Z] == distinct_types[j]) {
          log_lambda[i, j] <- 1
        }
      }
    }
  }

  rho_offset <- rep(0, n_Z + n_D, p)
  alpha_list <- vector(mode = "list", length = p * (p + 1) / 2)
  alpha_list <- lapply(alpha_list, function(x) rep(0, n_Z + n_D))

  # TODO: Horrible, vectorise
  for(i in seq_len(n_Z + n_D)) {
    if(i <= n_Z) {
      location <- c(configuration@x[i], configuration@y[i])
      type <- types(configuration)[i]
      type_index <- match(type, distinct_types)

      # TODO: We have i, should be quicker to remove from configuration
      disp <- compute_delta_phi_dispersion(remove_from_configuration(configuration, location, type), location, type_index - 1, p, model, radius)
    } else {
      location <- c(D@x[i - n_Z], D@y[i - n_Z])
      type_index <- match(types(D)[i - n_Z], distinct_types)

      disp <- compute_delta_phi_dispersion(configuration, location, type_index - 1, p, model, radius)
    }

    rho_offset[i] <- rho[type_index]

    index <- 1
    for(j in seq_len(p)) {
      for(k in j:p) {
        # TODO: Resetting names too many times in the `i` for loop.
        names(alpha_list)[index] <- paste0("alpha_", j, k)

        #alpha_list[[index]][i] <-  disp[j, k]
        if(j == type_index) {
          alpha_list[[index]][i] <- disp[k]
        } else if(k == type_index) {
          alpha_list[[index]][i] <- disp[j]
        }

        index <- index + 1
      }
    }
  }

  data <- as.data.frame(list(response = response, log_lambda = log_lambda, alpha_list, rho = rho_offset))

  formula <- paste("response ~ 0 + log_lambda + ", paste(names(alpha_list), collapse = " + "), " + offset(-log(rho))", sep = "")
  formula <- as.formula(formula)

  g <- glm(formula, data = data, family = binomial())

  if(print) {
    print(summary(g))
  }

  g
}