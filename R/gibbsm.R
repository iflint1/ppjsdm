#' Fit a multivariate Gibbs model to a dataset.
#'
#' @param configuration A configuration assumed to be a draw from the multivariate Gibbs.
#' @param window Observation window.
#' @param covariates Environmental covariates driving the intensity.
#' @param traits Species' traits.
#' @param model String to represent the model we're calibrating. You can check the currently authorised models with a call to `show_model()`.
#' @param radius Interaction radius.
#' @param print Print the fitting summary?
#' @importFrom stats as.formula binomial glm
#' @export
gibbsm <- function(configuration, window = Rectangle_window(), covariates = list(), traits = list(), model = "identity", radius = matrix(0), print = TRUE) {
  # num <- get_number_points(configuration)
  # num <- unlist(num, use.names = FALSE)
  # if(all(num == 0)) {
  #   stop("Empty configuration.")
  # }
  # # This is the guideline from the Baddeley et al. paper, see p. 8 therein.
  # rho_times_volume <- 4 * num
  # rho <- rho_times_volume / window_volume(window)
  #
  # types <- types(configuration)
  # distinct_types <- levels(types)
  # ntypes <- length(distinct_types)
  # ncovariates <- length(covariates)
  # ntraits <- length(traits)
  #
  # D <- rbinomialpp(window = window, n = rho_times_volume, types = distinct_types)
  #
  # n_Z <- get_number_points(configuration, total = TRUE)
  # n_D <- get_number_points(D, total = TRUE)
  # response <- c(rep.int(1, n_Z), rep.int(0, n_D))
  # log_lambda <- matrix(0, n_Z + n_D, ntypes)
  # rho_offset <- rep(0, n_Z + n_D, ntypes)
  #
  # alpha_list <- vector(mode = "list", length = ntypes * (ntypes + 1) / 2)
  # alpha_list <- lapply(alpha_list, function(x) rep.int(0, n_Z + n_D))
  #
  # covariate_list <- vector(mode = "list", length = ncovariates)
  # covariate_list <- lapply(covariate_list, function(x) rep.int(0, n_Z + n_D))
  # # TODO: When this is moved to C++, fill in blanks with type1 / type2, etc.
  # names(covariate_list) <- names(covariates)
  #
  # trait_list <- vector(mode = "list", length = ntraits)
  # trait_list <- lapply(trait_list, function(x) rep.int(0, n_Z + n_D))
  # # TODO: When this is moved to C++, fill in blanks with type1 / type2, etc.
  # names(trait_list) <- names(traits)
  #
  # index <- 1
  # for(j in seq_len(ntypes)) {
  #   for(k in j:ntypes) {
  #     names(alpha_list)[index] <- paste0("alpha", j, "_", k)
  #     index <- index + 1
  #   }
  # }
  #
  # # TODO: Might be able to vectorise -- Not urgent since this will be converted to C++ in the longer term.
  # for(i in seq_len(n_Z + n_D)) {
  #   if(i <= n_Z) {
  #     location <- c(x_coordinates(configuration)[i], y_coordinates(configuration)[i])
  #     type <- types(configuration)[i]
  #     type_index <- match(type, distinct_types)
  #
  #     dispersion <- compute_delta_phi_dispersion(remove_from_configuration_by_index(configuration, i), location, type_index - 1, ntypes, model, radius)
  #   } else {
  #     location <- c(x_coordinates(D)[i - n_Z], y_coordinates(D)[i - n_Z])
  #     type <- types(D)[i - n_Z]
  #     type_index <- match(type, distinct_types)
  #
  #     dispersion <- compute_delta_phi_dispersion(configuration, location, type_index - 1, ntypes, model, radius)
  #   }
  #
  #   rho_offset[i] <- rho[type_index]
  #   for(j in seq_len(ncovariates)) {
  #     covariate_list[[j]][i] <- covariates[[j]](location[1], location[2])
  #   }
  #   for(j in seq_len(ntraits)) {
  #     tr <- traits[[j]][type_index, 1:ntypes]
  #     trait_list[[j]][i] <- tr %*% dispersion
  #   }
  #
  #   index <- 1
  #   for(j in seq_len(ntypes)) {
  #     if(type == distinct_types[j]) {
  #       log_lambda[i, j] <- 1
  #     }
  #     for(k in j:ntypes) {
  #       if(j == type_index) {
  #         alpha_list[[index]][i] <- dispersion[k]
  #       } else if(k == type_index) {
  #         alpha_list[[index]][i] <- dispersion[j]
  #       }
  #
  #       index <- index + 1
  #     }
  #   }
  # }

  ret <- prepare_gibbsm_data(configuration, window, covariates, traits, model, radius)
  covariate_list <- ret$covariates
  trait_list <- ret$traits
  alpha_list <- ret$alpha
  response <- ret$response
  log_lambda <- ret$log_lambda
  rho <- ret$rho
  formula <- ret$formula

  data <- cbind(covariate_list, trait_list, alpha_list, response, log_lambda, rho)

  g <- glm(as.formula(formula), data = as.data.frame(data), family = binomial())

  if(print) {
    print(summary(g))
  }

  g
}
