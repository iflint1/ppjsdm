# dispersion <- function (configuration, t = levels(types(configuration))) {
#   X <- cbind(configuration@x, configuration@y)
#   species <- types(configuration)
#
#   phi <- function(x) {
#     R <- 0.1
#     x_larger_than_R <- x > R
#     x_smaller_than_R <- x <= R
#
#     x[x_larger_than_R] <- 0
#     x[x_smaller_than_R] <- 1
#     x
#   }
#
#   #phi <- function(x) {
#   #  x
#   #}
#
#   # dense Euclidean distance matrix
#   dd <- as.matrix(dist(X))
#   dists <- as.matrix(apply(dd, 2, phi))
#
#   total_species <- length(t)
#
#   # Overall abundance vector
#   M <- rep(0, total_species)
#   MM <- get_number_points(configuration)
#   for (i in seq_len(total_species)) {
#     type_index <- match(t[i], levels(species))
#     if(!is.na(type_index)) {
#       M[i] <- MM[type_index]
#     }
#   }
#
#   # Disperson matrix
#   # TODO: Make this independent of unit rectangle.
#   default_dispersion <- phi(sqrt(2))
#   D <- matrix(default_dispersion, total_species, total_species)
#   for (i in seq_len(total_species)) {
#     for (j in seq_len(total_species)) {
#       dist_sub <- dists[species == t[i], species == t[j]]
#
#       if(M[i] > 0 && M[j] > 0) {
#         if (i == j) {
#           # only upper-triangular elements for self-interaction
#           dist_sub <- dist_sub[upper.tri(dist_sub)]
#           if(M[i] > 1) {
#             denom <- M[i] * (M[i] - 1) / 2
#           } else {
#             denom <- 1
#             stopifnot(sum(dist_sub) == 0)
#           }
#         } else {
#           denom <- M[i] * M[j]
#         }
#
#         # denom <- 1
#
#         # sum distances between relevant pairs of points
#         D[i, j] <- sum(dist_sub) / denom
#       }
#     }
#   }
#   D
# }

papangelou <- function(configuration, x, type, alpha, lambda, types) {
  # type_index <- match(type, types)
  #
  # stopifnot(!is.na(type_index))
  #
  # number_types <- length(types)
  #
  # D <- compute_phi_dispersion(configuration, number_types)
  # D_plus <- compute_phi_dispersion(add_to_configuration(configuration, x, type), number_types)
  #
  # alpha_i <- alpha[type_index, ]
  #
  # D_i <- (D_plus - D)[type_index, ]
  #
  # lambda[type_index] * exp(-alpha_i %*% D_i)

  type_index <- match(type, types)

  stopifnot(!is.na(type_index))

  number_types <- length(types)

  delta_D <- compute_delta_phi_dispersion(configuration, x, type_index - 1, number_types)

  alpha_i <- alpha[type_index, ]

  lambda[type_index] * exp(alpha_i %*% delta_D)
}


#' Sample from a multivariate Gibbs point process.
#'
#' @param alpha alpha
#' @param lambda lambda
#' @param nsim Number of steps in the Markov Chain algorithm
#' @param prob Probability of choosing a birth.
#' @param verbose Print extra information
#' @param return_full_chain Return the full chain or only the last element?
#' @importFrom stats runif
#' @export
rgibbs <- function (alpha, lambda, nsim = 30000, prob = 0.5, verbose = TRUE, return_full_chain = FALSE) {
  acceptance <- matrix(NA, nsim, 2)
  results <- list()

  current <- rppp(Rectangle_window(), lambda)

  # Save full list of species before going further in the algorithm; note that some species may be lost because of removal of points.
  types <- levels(types(current))

  for (i in 1:nsim) {
    U <- runif(2)
    if(U[1] <= prob) {
      # Uniformly distributed type
      type <- sample(types, 1)
      # Uniformly distributed location
      location <- runif(2)

      # TODO: Implement other windows of volume != 1
      volume <- 1
      birth_ratio <- papangelou(current, location, type, alpha, lambda, types) * (1 - prob) * volume * length(types) / (prob * (1 + sum(get_number_points(current))))

      if(U[2] <= birth_ratio) {
        current <- add_to_configuration(current, location, type)
        acceptance[i, 1] <- 1
      } else {
        acceptance[i, 1] <- 0
      }
    } else {
      if(!is_empty(current)) {
        current_number <- sum(get_number_points(current))
        removed <- remove_random_point(current)
        volume <- 1

        death_ratio <- prob * current_number / (length(types) * (1 - prob) * volume * papangelou(removed$configuration, removed$location, removed$type, alpha, lambda, types))
        if(U[2] <= death_ratio) {
          current <- removed$configuration
          acceptance[i, 2] <- 1
        } else {
          acceptance[i, 2] <- 0
        }
      }
    }

    results[[i]] <- current

  }

  if (verbose) {
    birth_acceptance <- acceptance[, 1]
    birth_acceptance <- birth_acceptance[!is.na(birth_acceptance)]

    msg <- sprintf("%i proposed births, %i accepted births (%s%s)",
                   length(birth_acceptance),
                   sum(birth_acceptance),
                   round(mean(birth_acceptance) * 100),
                   "%")
    print(msg)

    death_acceptance <- acceptance[, 2]
    death_acceptance <- death_acceptance[!is.na(death_acceptance)]

    msg <- sprintf("%i proposed deaths, %i accepted deaths (%s%s)",
                   length(death_acceptance),
                   sum(death_acceptance),
                   round(mean(death_acceptance) * 100),
                   "%")
    print(msg)
  }

  if(return_full_chain) {
    results
  } else {
    results[[nsim]]
  }
}
