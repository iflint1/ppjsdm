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
#' @param how Used when `type1` and/or `type2` are vectors of types, in which case we need to either combine them by pairs
#' (this assumes that they have the same length) or consider all possible combinations of `type1` with `type2`.
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
#' fit <- ppjsdm::gibbsm(configuration, short_range = list(matrix(0.05, 2, 2), matrix(0.1, 2, 2)), medium_range = matrix(0.2, 2, 2), long_range = matrix(0.3, 2, 2))
#'
#' plot(ppjsdm::potentials(fit))
#' plot(ppjsdm::potentials(fit, type1 = "A", type2 = "B"))
#' plot(ppjsdm::potentials(fit, type1 = c("A", "B"), type2 = c("A", "B")))
#'
#' # User wants to give names to the potentials
#'
#' fit <- ppjsdm::gibbsm(configuration, short_range = list(`0.05m potential` = matrix(0.05, 2, 2), `0.1m potential` = matrix(0.1, 2, 2)), medium_range = matrix(0.2, 2, 2), long_range = matrix(0.3, 2, 2))
#'
#' plot(ppjsdm::potentials(fit))
#'
potentials <- function(fit,
                       type1 = 1,
                       type2 = type1,
                       how = c("all", "pairs")) {
  # Check fit type
  if(!is(fit, "gibbsm")) {
    stop(paste0("The fit object does not have the right type, its class is ", class(fit)))
  }

  # Interpret the how argument
  how <- match.arg(how)

  # Compute how many pairs of types we need to consider
  to_consider <- if(how == "all") {
    length(type1) * (length(type2) + 1) / 2
  } else if(how == "pairs") {
    if(length(type1) != length(type2)) {
      stop("Cannot use this option if the length of type1 and type2 are different.")
    }
    length(type1)
  } else {
    stop("Unrecognised option.")
  }

  # Convert to list
  type1 <- as.list(type1)
  type2 <- as.list(type2)

  # Initialise return variables
  short_range <- vector(mode = "list", length = to_consider)
  long_range <- vector(mode = "list", length = to_consider)
  short_potentials <- vector(mode = "list", length = to_consider)
  medium_potential <- vector(mode = "list", length = to_consider)
  overall <- vector(mode = "list", length = to_consider)
  t1 <- vector(mode = "list", length = to_consider)
  t2 <- vector(mode = "list", length = to_consider)

  k1 <- k2 <- 1

  for(k in seq_len(to_consider)) {
    # Current types
    t1[[k]] <- if(is.numeric(type1[[k1]])) {
      rownames(fit$coefficients$medium_range)[type1[[k1]]]
    } else {
      type1[[k1]]
    }
    t2[[k]] <- if(is.numeric(type2[[k2]])) {
      rownames(fit$coefficients$medium_range)[type2[[k2]]]
    } else {
      type2[[k2]]
    }

    # Set indices
    if(how == "all") {
      k2 <- k2 + 1
      if(k2 > length(type2)) {
        k1 <- k1 + 1
        k2 <- k1
      }
    } else {
      k1 <- k2 <- k + 1
    }
  }

  eval_k <- function(k) {
    # Extract variables relevant to the short-range potentials
    model <- fit$parameters$model
    short_range <- lapply(fit$coefficients$short_range, function(s) s[t1[[k]], t2[[k]]])
    alpha <- lapply(fit$coefficients$alpha, function(a) a[t1[[k]], t2[[k]]])

    # Construct the short-range potentials
    short_potentials <- setNames(lapply(seq_along(model), function(i) {
      mod <- model[[i]]
      if(mod == "exponential") {
        function(x) alpha[[i]] * exp(-log(2.) * x / short_range[[i]])
      } else if(mod == "square_exponential") {
        function(x) alpha[[i]] * exp(-log(2.) * x^2 / short_range[[i]]^2)
      } else if(mod == "bump") {
        function(x) alpha[[i]] * (1 - exp(-short_range[[i]] * log(2.) / x))
      } else if(mod == "square_bump") {
        function(x) alpha[[i]] * (1 - exp(-short_range[[i]] * short_range[[i]] * log(2.) / (x * x)))
      } else if(mod == "Geyer") {
        function(x) alpha[[i]] * ifelse(x <= short_range[[i]], 1., 0.)
      } else if(mod == "linear") {
        function(x) alpha[[i]] * pmax(0, 1. - x / short_range[[i]])
      } else {
        stop(paste0("Short-range model not recognised: ", mod))
      }
    }), nm = fit$potential_names)

    # Extract variables relevant to the medium-range potentials
    medium_range_model <- fit$parameters$medium_range_model
    medium_range <- fit$coefficients$medium_range[t1[[k]], t2[[k]]]
    long_range <- fit$coefficients$long_range[t1[[k]], t2[[k]]]
    gamma <- fit$coefficients$gamma[t1[[k]], t2[[k]]]

    # Construct the medium-range potentials
    medium_potential <- if(medium_range_model == "square_exponential") {
      function(x) gamma * exp(-4. * log(2.) * ((medium_range + long_range) / 2. - x)^2 /
                                (medium_range - long_range)^2)
    } else if(medium_range_model == "half_square_exponential") {
      function(x) gamma * ifelse(x > medium_range, exp(-log(2.) * (x - medium_range)^2 /
                                                         (long_range - medium_range)^2), 0.)
    } else if(medium_range_model == "Geyer") {
      function(x) gamma * ifelse(x <= long_range & x >= medium_range, 1., 0.)
    } else if(medium_range_model == "linear") {
      function(x) gamma * ifelse(2. * x <= medium_range + long_range,
                                 ifelse(x <= medium_range, 0., 2. /
                                          (long_range - medium_range) *
                                          (x - medium_range)),
                                 ifelse(x >= long_range, 0., 2. /
                                          (long_range - medium_range) *
                                          (long_range - x)))
    } else if(medium_range_model == "half_exponential") {
      function(x) gamma * ifelse(x >= medium_range, exp(-log(2.) * (x - medium_range) /
                                                          (long_range - medium_range)),
                                 0.)
    } else if(medium_range_model == "exponential") {
      function(x) gamma * exp(-2. * log(2.) * abs(x - 0.5 * (long_range + medium_range))
                              / (long_range - medium_range))
    } else if(medium_range_model == "bump") {
      function(x) {
        me <- medium_range
        hi <- long_range
        gamma * (1. - exp(-0.5 * sign(x - 0.5 * (me + hi)) * log(2.) * (hi - me) / (x - 0.5 * (me + hi))))
      }
    } else if(medium_range_model == "square_bump") {
      function(x) {
        me <- medium_range
        hi <- long_range
        gamma * (1.0 - exp(-0.25 * log(2.) * (hi - me)^2 / (x - 0.5 * (me + hi))^2))
      }
    } else if(medium_range_model == "tanh") {
      function(x) {
        me <- medium_range
        hi <- long_range
        gamma * (1. / (2. * tanh(5. / 2.)) * (tanh(5. / (hi - me) * (x - me)) + tanh(5. / (hi - me) * (hi - x))))
      }
    } else {
      stop(paste0("Medium-range model not recognised: ", medium_range_model))
    }

    # Compute overall potential
    overall <- function(x) Reduce("+", lapply(short_potentials, function(pot) pot(x))) + medium_potential(x)

    list(short = short_potentials,
         medium = medium_potential,
         overall = overall,
         short_range = short_range,
         long_range = long_range)
  }

  # Return object
  rev <- lapply(seq_len(to_consider), eval_k)
  ret <- list(short = lapply(rev, function(r) r$short),
              medium = lapply(rev, function(r) r$medium),
              overall = lapply(rev, function(r) r$overall),
              short_range = lapply(rev, function(r) r$short_range),
              long_range = lapply(rev, function(r) r$long_range),
              type1 = t1,
              type2 = t2)
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

#' @importFrom ggnewscale new_scale_colour
#' @importFrom ggplot2 aes element_blank geom_line ggplot guide_legend guides scale_colour_viridis_d scale_linetype_manual theme theme_minimal xlab ylab
#' @importFrom rlang .data
#' @importFrom stats complete.cases median
#' @method plot potentials
#' @export
plot.potentials <- function(x, base_size = 11, compute_average = TRUE, compute_median = TRUE, ...) {
  # Compute xseq
  max_distance <- 1.5 * max(max(sapply(x$short_range, function(y) Reduce("max", y))), Reduce("max", x$long_range))
  xseq <- seq(from = 0, to = max_distance, length.out = 1e3)

  # Number of pairs of types to consider
  npairs <- length(x$type1)

  df <- NULL

  for(k in seq_len(npairs)) {
    pair <- paste0(x$type1[[k]], enc2utf8(" \u2194 "), x$type2[[k]])
    df_current <- data.frame(x = xseq,
                             Overall = x$overall[[k]](xseq),
                             Medium = x$medium[[k]](xseq),
                             pair = pair)

    for(i in seq_len(length(x$short[[k]]))) {
      df_current[, paste0("short", i)] <- x$short[[k]][[i]](xseq)
      if(is.null(names(x$short[[k]])[i])) {
        df_current[, paste0("name", i)] <- paste0("Short ", i)
      } else if(names(x$short[[k]])[i] == "") { # Cannot put this as a conditional in the previous case due to no short-circuiting
        df_current[, paste0("name", i)] <- paste0("Short ", i)
      } else {
        df_current[, paste0("name", i)] <- names(x$short[[k]])[i]
      }
    }

    df <- rbind(df, df_current)
  }

  df <- df[complete.cases(df), ]

  is_one_pair <- length(unique(df$pair)) == 1

  if(!is_one_pair & (compute_average | compute_median)) {
    pairs <- unique(df$pair)

    avg <- lapply(pairs, function(pr) {
      df$Overall[df$pair == pr]
    })

    if(compute_median) {
      md <- sapply(seq_along(avg[[1]]), function(i) median(sapply(avg, function(a) a[i])))
    }

    if(compute_average) {
      avg <- Reduce("+", avg) / length(avg)
    }
  }

  g <- ggplot(data = df) +
    geom_line(aes(x = .data$x, y = 0), colour = "black") +
    xlab(NULL) +
    ylab(NULL) +
    theme_minimal(base_size = base_size) +
    theme(legend.title = element_blank())

  if(is_one_pair) {
    g <- g + geom_line(aes(x = .data$x, y = .data$Overall, colour = "Overall"), size = 1.5)
  } else {
    g <- g + geom_line(aes(x = .data$x, y = .data$Overall, linetype = "Overall", colour = .data$pair), size = 1.5, alpha = 0.6)
  }


  nshorts <- sum(startsWith(colnames(df), "short"))
  for(i in seq_len(nshorts)) {
    if(is_one_pair) {
      g <- g + geom_line(aes(x = .data$x, y = .data[[paste0("short", i)]],
                             colour = .data[[paste0("name", i)]]), size = 1, alpha = 0.7)
    } else {
      g <- g + geom_line(aes(x = .data$x, y = .data[[paste0("short", i)]],
                             linetype = "Short/Medium", colour = .data$pair), size = 1, alpha = 0.5)
    }
  }

  if(any(df$Medium != 0, na.rm = TRUE)) {
    if(is_one_pair) {
      g <- g + geom_line(aes(x = .data$x, y = .data$Medium, colour = "Medium"), size = 1, alpha = 0.7)
    } else {
      g <- g + geom_line(aes(x = .data$x, y = .data$Medium, linetype = "Short/Medium", colour = .data$pair), size = 1, alpha = 0.5)
    }
  }

  g <- g +
    scale_linetype_manual(values = c(Overall = "solid", `Short/Medium` = "dotted")) +
    scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "turbo",
                           guide = guide_legend(order = 2, nrow = 10)) +
    guides(linetype = guide_legend(order = 3,
                                   keywidth = 3,
                                   override.aes = list(linewidth = 1))) +
    new_scale_colour()

  if(!is_one_pair & compute_median) {
    g <- g + geom_line(data = data.frame(x = xseq, md = md), aes(x = .data$x, y = .data$md, colour = "Median"), alpha = 0.9, size = 2.5)
  }
  if(!is_one_pair & compute_average) {
    g <- g + geom_line(data = data.frame(x = xseq, avg = avg), aes(x = .data$x, y = .data$avg, colour = "Average"), alpha = 0.9, size = 2.5)
  }

  g <- g + scale_colour_manual(values = c(Average = "black", Median = "red"), guide = guide_legend(order = 1))

  g
}
