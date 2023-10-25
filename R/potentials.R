#' Recover the interaction potentials between two given types from a fitted `gibbsm` object.
#'
#' After calling gibbsm, this function allows the user to visualise the fitted
#' potentials. The function returns an object that can then be printed or plotted.
#'
#' @param fit A fit object returned by a call to `gibbsm`.
#' @param type1 First interacting types. Can be either an integer representing the type index,
#' or a string representing its label.
#' @param type2 Second interacting types. Can be either an integer representing the type index,
#' or a string representing its label.
#' @export
#' @md
#' @examples
#' set.seed(1)
#'
#' # Simulate a Poisson point process with two types
#'
#' configuration <- ppjsdm::rppp(lambda = c(A = 100, B = 100))
#'
#' # Fit with gibbsm
#'
#' fit <- ppjsdm::gibbsm(configuration, short_range = c(matrix(0.05, 2, 2), matrix(0.1, 2, 2)), medium_range = matrix(0.2, 2, 2), long_range = matrix(0.3, 2, 2))
#'
#' plot(ppjsdm::potentials(fit))
#' plot(ppjsdm::potentials(fit, type1 = "A", type2 = "B"))
#'
potentials <- function(fit,
                       type1 = 1,
                       type2 = type1) {
  # Check fit type
  if(!is(fit, "gibbsm")) {
    stop(paste0("The fit object does not have the right type, its class is ", class(fit)))
  }

  # Extract variables relevant to the short-range potentials
  model <- fit$parameters$model
  short_range <- fit$coefficients$short_range
  alpha <- fit$coefficients$alpha

  # Construct the short-range potentials
  short_potentials <- lapply(seq_len(length(model)), function(i) {
    mod <- model[[i]]
    if(mod == "exponential") {
      function(x) alpha[[i]][type1, type2] * exp(-log(2) * x / short_range[[i]][type1, type2])
    } else if(mod == "square_exponential") {
      function(x) alpha[[i]][type1, type2] * exp(-log(2) * x^2 / short_range[[i]][type1, type2]^2)
    } else if(mod == "bump") {
      function(x) alpha[[i]][type1, type2] * (1 - exp(-short_range[[i]][type1, type2] * log(2) / x))
    } else if(mod == "square_bump") {
      function(x) alpha[[i]][type1, type2] * (1 - exp(-short_range[[i]][type1, type2] *
                                                              short_range[[i]][type1, type2] * log(2) / (x * x)))
    } else if(mod == "Geyer") {
      function(x) alpha[[i]][type1, type2] * ifelse(x <= short_range[[i]][type1, type2], 1, 0)
    } else if(mod == "linear") {
      function(x) alpha[[i]][type1, type2] * pmax(0, 1. - x / short_range[[i]][type1, type2])
    } else {
      stop(paste0("Short-range model not recognised: ", mod))
    }
  })

  # Extract variables relevant to the medium-range potentials
  medium_range_model <- fit$parameters$medium_range_model
  medium_range <- fit$coefficients$medium_range[type1, type2]
  long_range <- fit$coefficients$long_range[type1, type2]
  gamma <- fit$coefficients$gamma[type1, type2]

  # Construct the medium-range potentials
  medium_potential <- if(medium_range_model == "square_exponential") {
    function(x) gamma * exp(-4 * log(2) * ((medium_range + long_range) / 2 - x)^2 /
                      (medium_range - long_range)^2)
  } else if(medium_range_model == "half_square_exponential") {
    function(x) gamma * ifelse(x > medium_range, exp(-log(2) * (x - medium_range)^2 /
                                                                   (long_range - medium_range)^2), 0.)
  } else if(medium_range_model == "Geyer") {
    function(x) gamma * ifelse(x <= long_range & x >= medium_range, 1, 0)
  } else if(medium_range_model == "linear") {
    function(x) gamma * ifelse(2 * x <= medium_range + long_range,
                       ifelse(x <= medium_range, 0., 2. /
                                (long_range - medium_range) *
                                (x - medium_range)),
                       ifelse(x >= long_range, 0., 2. /
                                (long_range - medium_range) *
                                (long_range - x)))
  } else if(medium_range_model == "half_exponential") {
    function(x) gamma * ifelse(x >= medium_range, exp(-log(2) *
                                                                    (x - medium_range) /
                                                                    (long_range - medium_range)),
                       0.)
  } else if(medium_range_model == "exponential") {
    function(x) gamma * exp(-2 * log(2) * abs(x - 0.5 * (long_range + medium_range))
                    / (long_range - medium_range))
  } else if(medium_range_model == "bump") {
    function(x) {
      me <- medium_range
      hi <- long_range
      gamma * (1.0 - exp(-0.5 * sign(x - 0.5 * (me + hi)) * log(2) * (hi - me) / (x - 0.5 * (me + hi))))
    }
  } else if(medium_range_model == "square_bump") {
    function(x) {
      me <- medium_range
      hi <- long_range
      gamma * (1.0 - exp(-0.25 * log(2) * (hi - me)^2 / (x - 0.5 * (me + hi))^2))
    }
  } else if(medium_range_model == "tanh") {
    function(x) {
      me <- medium_range
      hi <- long_range
      gamma * (1 / (2 * tanh(5 / 2)) * (tanh(5 / (hi - me) * (x - me)) + tanh(5 / (hi - me) * (hi - x))))
    }
  } else {
    stop(paste0("Medium-range model not recognised: ", medium_range_model))
  }

  # Return object
  ret <- list(short = short_potentials,
              medium = medium_potential,
              overall = function(x) Reduce("+", lapply(short_potentials, function(pot) pot(x))) +
                medium_potential(x),
              long_range = long_range,
              short_range = sapply(seq_len(length(model)), function(i) short_range[[i]][type1, type2]))
  class(ret) <- "potentials"
  ret
}

format_potentials <- function(potentials) {
  str <- paste0("Interaction potentials corresponding to a provided fit.\n\n",
                "The object contains three slots:\n- obj$short is a list containing functions of x, ",
                "each entry corresponding to one of the short-range potentials.\n",
                "- obj$medium contains a function of x representing the medium-range potential.\n",
                "- obj$overall contains a function of x representing the sum of all potentials.")
  str
}

#' @method print potentials
#' @export
print.potentials <- function(x, ...) {
  cat(format_potentials(x))
}

#' @importFrom ggplot2 aes aes_string element_blank geom_line ggplot theme theme_minimal xlab ylab
#' @importFrom stats complete.cases
#' @method plot potentials
#' @export
plot.potentials <- function(x, ...) {
  max_distance <- 2 * max(Reduce("max", x$short_range), x$long_range)
  xseq <- seq(from = 0, to = max_distance, length.out = 1e3)
  df <- data.frame(x = xseq,
                   overall = x$overall(xseq),
                   medium = x$medium(xseq))

  for(i in seq_len(length(x$short))) {
    df[, paste0("short", i)] <- x$short[[i]](xseq)
  }

  df <- df[complete.cases(df), ]

  Overall <- "Overall"
  g <- ggplot(data = df) +
    geom_line(aes_string(x = "x", y = "overall", colour = "Overall"), size = 2) +
    geom_line(aes(x = x, y = 0), colour = "black") +
    xlab(NULL) +
    ylab(NULL) +
    theme_minimal() +
    theme(legend.title = element_blank())

  for(i in seq_len(length(x$short))) {
    assign(paste0("name", i), paste0("Short ", i))
    g <- g + geom_line(aes_string(x = "x", y = paste0("short", i), colour = paste0("name", i)), size = 1, alpha = 0.8)
  }

  if(any(df$medium != 0, na.rm = TRUE)) {
    Medium <- "Medium"
    g <- g + geom_line(aes_string(x = "x", y = "medium", colour = "Medium"), size = 1, alpha = 0.8)
  }
  g
}
