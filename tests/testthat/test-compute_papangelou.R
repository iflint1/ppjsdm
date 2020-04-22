library(ppjsdm)
library(spatstat)

context("compute_papangelou")

test_that("Correct Papangelou conditional intensity value", {
  set.seed(42)
  Ntests <- 100
  for(nt in seq_len(Ntests)) {
    if(nt %% 3 == 0) {
      varphi <- function(x, i, j) ifelse(x <= r_1[i, j], 1, 0)
      model <- "Geyer"
    } else if(nt %% 3 == 1) {
      varphi <- function(x, i, j) exp(-log(2) * x / r_1[i, j])
      model <- "exponential"
    } else if(nt %% 3 == 2) {
      varphi <- function(x, i, j) 1 - exp(-r_1[i, j] * r_1[i, j] * log(2) / (x * x))
      model <- "square_bump"
    }

    if(nt %% 2 == 0) {
      psi <- function(x, i, j) ifelse(x <= r_3[i, j] & x >= r_2[i, j], 1, 0)
      medium_range_model <- "Geyer"
    } else if(nt %% 2 == 1) {
      psi <- function(x, i, j) exp(-4 * log(2) * ((r_2[i, j] + r_3[i, j]) / 2 - x)^2 / (r_2[i, j] - r_3[i, j])^2)
      medium_range_model <- "square_exponential"
    }

    configuration <- Configuration(x = runif(20, 0, 1), y = runif(20, 0, 1), types = c(rep("a", 9), rep("b", 11)))
    configuration_a <- Configuration(x = x_coordinates(configuration)[types(configuration) == "a"], y = y_coordinates(configuration)[types(configuration) == "a"], types = types(configuration)[types(configuration) == "a"])
    configuration_b <- Configuration(x = x_coordinates(configuration)[types(configuration) == "b"], y = y_coordinates(configuration)[types(configuration) == "b"], types = types(configuration)[types(configuration) == "b"])
    lambda <- c(1, 2)
    alpha <- cbind(c(0.5, 1), c(1, -0.5))
    gamma <- cbind(c(2, 0), c(0, -2))
    r_1 <- matrix(runif(4, 0, 1), 2, 2)
    r_1 <- r_1 %*% t(r_1)
    r_2 <- matrix(runif(4, 0, 1), 2, 2)
    r_2 <- r_2 %*% t(r_2)
    r_3 <- matrix(runif(4, 0, 1), 2, 2)
    r_3 <- r_3 %*% t(r_3)

    beta <- matrix(4, 2, 1)
    covariate <- function(x, y) x + y
    x <- runif(2, 0, 1)
    window <- owin(c(x[1], x[1]), c(x[2], x[2]))
    # IMPORTANT NOTE: I'm assuming in places below that type = 1, so you can't change type below and expect everything to work.
    type <- 1

    N <- 2
    u_term_fun <- function(x, conf, i, j) {
      s <- sapply(seq_len(length(conf$x)), function(k) varphi(dist(rbind(x, c(conf$x[k], conf$y[k]))), i, j))
      sum(s[order(s, decreasing = TRUE)[1:N]])
    }
    v_term_fun <- function(x, conf, i, j) {
      s <- sapply(seq_len(length(conf$x)), function(k) psi(dist(rbind(x, c(conf$x[k], conf$y[k]))), i, j))
      sum(s[order(s, decreasing = TRUE)[1:N]])
    }

    u_term1 <- u_term_fun(x, configuration_a, type, 1) + sum(sapply(seq_len(length(configuration_a$x)), function(i) {
      configuration_minus <- Configuration(x = configuration_a$x[-i], y = configuration_a$y[-i], types = configuration_a$types[-i])
      configuration_minus_plus_x <- Configuration(x = c(configuration_a$x[-i], x[1]), y = c(configuration_a$y[-i], x[2]), types = configuration_a$types)
      p <- c(configuration_a$x[i], configuration_a$y[i])
      z <- u_term_fun(p, configuration_minus_plus_x, 1, type) - u_term_fun(p, configuration_minus, 1, type)
      z
    }))
    u_term2 <- u_term_fun(x, configuration_b, type, 2) + sum(sapply(seq_len(length(configuration_b$x)), function(i) {
      configuration_minus <- configuration_a
      configuration_minus_plus_x <- Configuration(x = c(configuration_a$x, x[1]), y = c(configuration_a$y, x[2]), types = factor(c(configuration_a$types, 1), labels = "a"))
      p <- c(configuration_b$x[i], configuration_b$y[i])
      z <- u_term_fun(p, configuration_minus_plus_x, 2, type) - u_term_fun(p, configuration_minus, 2, type)
      z
    }))
    u_term <- c(u_term1, u_term2)

    v_term1 <- v_term_fun(x, configuration_a, type, 1) + sum(sapply(seq_len(length(configuration_a$x)), function(i) {
      configuration_minus <- Configuration(x = configuration_a$x[-i], y = configuration_a$y[-i], types = configuration_a$types[-i])
      configuration_minus_plus_x <- Configuration(x = c(configuration_a$x[-i], x[1]), y = c(configuration_a$y[-i], x[2]), types = configuration_a$types)
      p <- c(configuration_a$x[i], configuration_a$y[i])
      z <- v_term_fun(p, configuration_minus_plus_x, 1, type) - v_term_fun(p, configuration_minus, 1, type)
      z
    }))
    v_term2 <- v_term_fun(x, configuration_b, type, 2) + sum(sapply(seq_len(length(configuration_b$x)), function(i) {
      configuration_minus <- configuration_a
      configuration_minus_plus_x <- Configuration(x = c(configuration_a$x, x[1]), y = c(configuration_a$y, x[2]), types = factor(c(configuration_a$types, 1), labels = "a"))
      p <- c(configuration_b$x[i], configuration_b$y[i])
      z <- v_term_fun(p, configuration_minus_plus_x, 2, type) - v_term_fun(p, configuration_minus, 2, type)
      z
    }))
    v_term <- c(v_term1, v_term2)

    papangelou_direct <- as.double(lambda[type] * exp(beta[type] * covariate(x[1], x[2]) + alpha[type, ] %*% u_term + gamma[type, ] %*% v_term))
    papangelou_package <- ppjsdm::compute_papangelou(configuration = configuration,
                                                     x = x[1],
                                                     y = x[2],
                                                     type = type,
                                                     model = model,
                                                     medium_range_model = medium_range_model,
                                                     alpha = alpha,
                                                     lambda = lambda,
                                                     beta = beta,
                                                     gamma = gamma,
                                                     covariates = list(cov = as.im(covariate, W = window)),
                                                     short_range = r_1,
                                                     medium_range = r_2,
                                                     long_range = r_3,
                                                     saturation = N)
    expect_equal(papangelou_direct, papangelou_package)
  }
})
