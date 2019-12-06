p <- 2

Z <- get_csiro_trees_configuration(6, prevalent = p)

plot(Z)

# Z <- rgibbs(cbind(c(1, 1), c(1, 1)), c(10, 10))

window <- Rectangle_window(c(0, 100), c(0, 50))
radius <- 5

# This is the guideline from the Baddeley et al. paper
rho <- 4 * get_number_points(Z) / window_volume(window)

types <- levels(types(Z))

D <- rppp(window, lambda = rho, types = types)

n_Z <- sum(get_number_points(Z))
n_D <- sum(get_number_points(D))
response <- c(rep.int(1, n_Z), rep.int(0, n_D))

log_lambda <- matrix(0, n_Z + n_D, p)
# TODO: Horrible, vectorise
for(i in seq_len(n_Z + n_D)) {
  for(j in seq_len(p)) {
    if(i <= n_Z) {
      if(types(Z)[i] == types[j]) {
        log_lambda[i, j] <- 1
      }
    } else {
      if(types(D)[i - n_Z] == types[j]) {
        log_lambda[i, j] <- 1
      }
    }
  }
}

rho_offset <- rep(0, n_Z + n_D, p)
# TODO: Horrible, vectorise
for(i in seq_len(n_Z + n_D)) {
  if(i <= n_Z) {
    type <- types(Z)[i]
  } else {
    type <- types(D)[i - n_Z]
  }
  rho_offset[i] <- rho[match(type, types)]
}

alpha_list <- vector(mode = "list", length = p * (p + 1) / 2)
alpha_list <- lapply(alpha_list, function(x) rep(0, n_Z + n_D))

# alpha <- matrix(NA, n_Z + n_D, p * (p + 1) / 2)

# TODO: Horrible, vectorise
for(i in seq_len(n_Z + n_D)) {
  if(i <= n_Z) {
    location <- c(Z@x[i], Z@y[i])
    type <- types(Z)[i]

    # TODO: We have i, should be quicker to remove from configuration
    disp <- compute_phi_dispersion(remove_from_configuration(Z, location, type), p, radius = radius) - compute_phi_dispersion(Z, p, radius = radius)
  } else {
    location <- c(D@x[i - n_Z], D@y[i - n_Z])
    type <- types(D)[i - n_Z]

    disp <- compute_phi_dispersion(Z, p, radius = radius) - compute_phi_dispersion(add_to_configuration(Z, location, type), p, radius = radius)
  }

  index <- 1
  for(j in seq_len(p)) {
    for(k in j:p) {
      # TODO: Resetting names too many times
      names(alpha_list)[index] <- paste0("alpha_", j, k)
      alpha_list[[index]][i] <-  disp[j, k]
      index <- index + 1
    }
  }
  # index <- 1
  # for(j in seq_len(p)) {
  #   for(k in j:p) {
  #     alpha[i, index] <- disp[j, k]
  #     index <- index + 1
  #   }
  # }
}

data <- as.data.frame(list(response = response, log_lambda = log_lambda, alpha_list, rho = rho_offset))

formula <- paste("response ~ 0 + log_lambda + ", paste(names(alpha_list), collapse = " + "), " + offset(-log(rho))", sep = "")
formula <- as.formula(formula)

g <- glm(formula, data = data, family = binomial())

print(summary(g))
