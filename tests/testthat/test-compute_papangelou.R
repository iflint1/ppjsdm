library(ppjsdm)

context("compute_papangelou")

test_that("Default values", {
  set.seed(1)
  x <- stats::runif(n = 1)
  y <- stats::runif(n = 1)
  xs <- stats::runif(n = 2)
  ys <- stats::runif(n = 2)
  marks <- stats::runif(n = 2)
  configuration <- ppjsdm::Configuration(x = xs, y = ys, types = c("a", "b"), marks = marks)
  papangelou <- ppjsdm::compute_papangelou(configuration = configuration,
                                           x = x,
                                           y = y)
  papangelou_explicit <- ppjsdm::compute_papangelou(configuration = configuration,
                                                    x = x,
                                                    y = y,
                                                    type = 1,
                                                    mark = 1.0,
                                                    model = "square_bump",
                                                    medium_range_model = "square_exponential",
                                                    alpha = matrix(0, ncol = 2, nrow = 2),
                                                    beta0 = c(0, 0),
                                                    beta = matrix(0, ncol = 0, nrow = 2),
                                                    gamma = matrix(0, ncol = 2, nrow = 2),
                                                    covariates = list(),
                                                    short_range = matrix(0, ncol = 2, nrow = 2),
                                                    medium_range = matrix(0, ncol = 2, nrow = 2),
                                                    long_range = matrix(0, ncol = 2, nrow = 2),
                                                    nthreads = 4,
                                                    saturation = 2)
  expect_equal(papangelou, papangelou_explicit)
})

test_that("Correct Papangelou conditional intensity value", {
  set.seed(42)
  Ntests <- 100
  for(nt in seq_len(Ntests)) {
    possible_short <- c("exponential",
                        "square_exponential",
                        "bump",
                        "square_bump",
                        "Geyer",
                        "linear")
    model <- sample(possible_short, 1)
    if(model == possible_short[1]) {
      varphi <- function(x, i, j) exp(-log(2) * x / r_1[i, j])
    } else if(model == possible_short[2]) {
      varphi <- function(x, i, j) exp(-log(2) * x^2 / r_1[i, j]^2)
    } else if(model == possible_short[3]) {
      varphi <- function(x, i, j) 1 - exp(-r_1[i, j] * log(2) / x)
    } else if(model == possible_short[4]) {
      varphi <- function(x, i, j) 1 - exp(-r_1[i, j] * r_1[i, j] * log(2) / (x * x))
    } else if(model == possible_short[5]) {
      varphi <- function(x, i, j) ifelse(x <= r_1[i, j], 1, 0)
    } else if(model == possible_short[6]) {
      varphi <- function(x, i, j) pmax(0, 1. - x / r_1[i, j])
    }

    possible_medium <- c("square_exponential",
                         "half_square_exponential",
                         "Geyer",
                         "linear",
                         "half_exponential",
                         "exponential",
                         "bump",
                         "square_bump",
                         "tanh")
    medium_range_model <- sample(possible_medium, 1)
    if(medium_range_model == possible_medium[1]) {
      psi <- function(x, i, j) exp(-4 * log(2) * ((r_2[i, j] + r_3[i, j]) / 2 - x)^2 / (r_2[i, j] - r_3[i, j])^2)
    } else if(medium_range_model == possible_medium[2]) {
      psi <- function(x, i, j) ifelse(x > r_2[i, j], exp(-log(2) * (x - r_2[i, j])^2 / (r_3[i, j] - r_2[i, j])^2), 0.)
    } else if(medium_range_model == possible_medium[3]) {
      psi <- function(x, i, j) ifelse(x <= r_3[i, j] & x >= r_2[i, j], 1, 0)
    } else if(medium_range_model == possible_medium[4]) {
      psi <- function(x, i, j) ifelse(2 * x <= r_2[i, j] + r_3[i, j],
                                      ifelse(x <= r_2[i, j], 0., 2. / (r_3[i, j] - r_2[i, j]) * (x - r_2[i, j])),
                                      ifelse(x >= r_3[i, j], 0., 2. / (r_3[i, j] - r_2[i, j]) * (r_3[i, j] - x)))
    } else if(medium_range_model == possible_medium[5]) {
      psi <- function(x, i, j) ifelse(x >= r_2[i, j], exp(-log(2) * (x - r_2[i, j]) / (r_3[i, j] - r_2[i, j])), 0.)
    } else if(medium_range_model == possible_medium[6]) {
      psi <- function(x, i, j) exp(-2 * log(2) * abs(x - 0.5 * (r_3[i, j] + r_2[i, j])) / (r_3[i, j] - r_2[i, j]))
    } else if(medium_range_model == possible_medium[7]) {
      psi <- function(x, i, j) {
        me <- r_2[i, j]
        hi <- r_3[i, j]
        1.0 - exp(-0.5 * sign(x - 0.5 * (me + hi)) * log(2) * (hi - me) / (x - 0.5 * (me + hi)))
      }
    } else if(medium_range_model == possible_medium[8]) {
      psi <- function(x, i, j) {
        me <- r_2[i, j]
        hi <- r_3[i, j]
        1.0 - exp(-0.25 * log(2) * (hi - me)^2 / (x - 0.5 * (me + hi))^2)
      }
    } else if(medium_range_model == possible_medium[9]) {
      psi <- function(x, i, j) {
        me <- r_2[i, j]
        hi <- r_3[i, j]
        1 / (2 * tanh(5 / 2)) * (tanh(5 / (hi - me) * (x - me)) + tanh(5 / (hi - me) * (hi - x)))
      }
    }

    configuration <- ppjsdm::Configuration(x = runif(30, 0, 1), y = runif(30, 0, 1), types = c(rep("a", 14), rep("b", 16)))
    configuration_by_type <- lapply(as.character(levels(ppjsdm::types(configuration))), function(type) {
      ppjsdm::Configuration(x = ppjsdm::x_coordinates(configuration)[ppjsdm::types(configuration) == type],
                            y = ppjsdm::y_coordinates(configuration)[ppjsdm::types(configuration) == type],
                            types = ppjsdm::types(configuration)[ppjsdm::types(configuration) == type])
    })

    beta0 <- c(1, 2)
    alpha <- cbind(c(0.5, 1), c(1, -0.5))
    gamma <- cbind(c(2, 0), c(0, -2))
    r_1 <- matrix(runif(4, 0, 0.1), 2, 2)
    r_1 <- r_1 %*% t(r_1)
    r_2 <- matrix(runif(4, 0, 0.1), 2, 2)
    r_2 <- r_2 %*% t(r_2)
    r_3 <- matrix(runif(4, 0, 0.1), 2, 2)
    r_3 <- r_2 + r_3 %*% t(r_3)

    N <- 2
    beta <- matrix(4, 2, 1)
    covariate <- function(x, y) x + y
    if(sample(c(TRUE, FALSE), 1)) { # Compute the Papangelou conditional intensity on a point of the configuration
      x <- c(configuration$x[1], configuration$y[1])
    } else { # Compute the Papangelou conditional intensity on an outside point
      x <- stats::runif(2, 0, 1)
    }
    window <- spatstat.geom::owin(c(0, 1), c(0, 1))
    type <- 1
    not_type <- 2

    # If x is in the configuration, remove it.
    x_is_in_configuration <- which(x[1] == configuration_by_type[[type]]$x & x[2] == configuration_by_type[[type]]$y)
    if(length(x_is_in_configuration) > 0) {
      configuration_by_type[[type]] <- ppjsdm::Configuration(x = configuration_by_type[[type]]$x[-x_is_in_configuration],
                                                             y = configuration_by_type[[type]]$y[-x_is_in_configuration],
                                                             types = configuration_by_type[[type]]$types[-x_is_in_configuration])
    }

    dispersion_term <- function(x, conf, i, j, range) {
      s <- sapply(seq_len(length(conf$x)), function(k) {
        if(range == "short") {
          varphi(dist(rbind(x, c(conf$x[k], conf$y[k]))), i, j)
        } else if(range == "medium") {
          psi(dist(rbind(x, c(conf$x[k], conf$y[k]))), i, j)
        } else {
          stop("Incorrect range value in dispersion term.")
        }
      })
      sum(s[order(s, decreasing = TRUE)[1:N]])
    }

    compute_dispersion_term <- function(range) {
      terms <- vector(mode = 'list', length = nlevels(ppjsdm::types(configuration)))
      terms[[type]] <- dispersion_term(x, configuration_by_type[[type]], type, type, range = range) +
        sum(sapply(seq_len(length(configuration_by_type[[type]]$x)), function(i) {
          configuration_minus <- ppjsdm::Configuration(x = configuration_by_type[[type]]$x[-i],
                                                       y = configuration_by_type[[type]]$y[-i],
                                                       types = configuration_by_type[[type]]$types[-i])
          configuration_minus_plus_x <- ppjsdm::Configuration(x = c(configuration_by_type[[type]]$x[-i], x[1]),
                                                              y = c(configuration_by_type[[type]]$y[-i], x[2]),
                                                              types = factor(c(configuration_by_type[[type]]$types[-i], type), labels = "a"))
          p <- c(configuration_by_type[[type]]$x[i], configuration_by_type[[type]]$y[i])
          dispersion_term(p, configuration_minus_plus_x, type, type, range = range) - dispersion_term(p, configuration_minus, type, type, range = range)
        }))

      terms[[not_type]] <- dispersion_term(x, configuration_by_type[[not_type]], type, not_type, range = range) + sum(sapply(seq_len(length(configuration_by_type[[not_type]]$x)), function(i) {
        configuration_minus <- configuration_by_type[[type]]
        configuration_minus_plus_x <- ppjsdm::Configuration(x = c(configuration_by_type[[type]]$x, x[1]),
                                                            y = c(configuration_by_type[[type]]$y, x[2]),
                                                            types = factor(c(configuration_by_type[[type]]$types, type), labels = "a"))
        p <- c(configuration_by_type[[not_type]]$x[i], configuration_by_type[[not_type]]$y[i])
        dispersion_term(p, configuration_minus_plus_x, not_type, type, range = range) - dispersion_term(p, configuration_minus, not_type, type, range = range)
      }))
      Reduce(c, terms)
    }

    u_term <- compute_dispersion_term(range = "short")
    v_term <- compute_dispersion_term(range = "medium")

    papangelou_direct <- as.double(exp(beta0[type] + beta[type] * covariate(x[1], x[2]) + alpha[type, ] %*% u_term + gamma[type, ] %*% v_term))
    papangelou_package <- ppjsdm::compute_papangelou(configuration = configuration,
                                                     x = x[1],
                                                     y = x[2],
                                                     type = type,
                                                     model = model,
                                                     medium_range_model = medium_range_model,
                                                     alpha = alpha,
                                                     beta0 = beta0,
                                                     beta = beta,
                                                     gamma = gamma,
                                                     covariates = list(cov = spatstat.geom::as.im(covariate(x[1], x[2]), W = window)),
                                                     short_range = r_1,
                                                     medium_range = r_2,
                                                     long_range = r_3,
                                                     nthreads = 4,
                                                     saturation = N)
    expect_equal(papangelou_direct, papangelou_package)
  }
})
