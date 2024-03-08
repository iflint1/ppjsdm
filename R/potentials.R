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
#' @param exclude_same_types Should pairs of potentials involving the same types be excluded from consideration?
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
                       how = c("all", "pairs"),
                       exclude_same_types = FALSE) {
  # Check fit type
  if(!is(fit, "gibbsm")) {
    stop(paste0("The fit object does not have the right type, its class is ", class(fit)))
  }

  # Interpret the how argument
  how <- match.arg(how)

  # Compute how many pairs of types we need to consider
  to_consider <- if(how == "all") {
    length(type1) * length(type2)
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

  t1 <- vector(mode = "character", length = to_consider)
  t2 <- vector(mode = "character", length = to_consider)

  k1 <- k2 <- 1

  for(k in seq_len(to_consider)) {
    # Current types
    t1[k] <- if(is.numeric(type1[[k1]])) {
      rownames(fit$coefficients$medium_range)[type1[[k1]]]
    } else {
      type1[[k1]]
    }
    t2[k] <- if(is.numeric(type2[[k2]])) {
      rownames(fit$coefficients$medium_range)[type2[[k2]]]
    } else {
      type2[[k2]]
    }

    # Set indices
    if(how == "all") {
      k2 <- k2 + 1
      if(k2 > length(type2)) {
        k1 <- k1 + 1
        k2 <- 1
      }
    } else {
      k1 <- k2 <- k + 1
    }
  }

  # Order the types by alphabetical order
  m <- pmin(t1, t2)
  M <- pmax(t1, t2)
  t1 <- m # Need to do this in two steps otherwise t1 is overwritten
  t2 <- M
  prs <- paste0(t1, t2)
  t1 <- t1[match(unique(prs), prs)]
  t2 <- t2[match(unique(prs), prs)]

  # Take care of exclude_same_types
  if(exclude_same_types) {
    z <- t1[t1 != t2] # Need to do this in two steps otherwise t1 is overwritten
    t2 <- t2[t1 != t2]
    t1 <- z
  }

  # Initialise return variables
  short_range <- vector(mode = "list", length = length(t1))
  long_range <- vector(mode = "list", length = length(t1))
  short_potentials <- vector(mode = "list", length = length(t1))
  medium_potential <- vector(mode = "list", length = length(t1))
  overall <- vector(mode = "list", length = length(t1))

  eval_k <- function(k) {
    # Extract variables relevant to the short-range potentials
    model <- fit$parameters$model
    short_range <- lapply(fit$coefficients$short_range, function(s) s[t1[k], t2[k]])
    alpha <- lapply(fit$coefficients$alpha, function(a) a[t1[k], t2[k]])

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
    medium_range <- fit$coefficients$medium_range[t1[k], t2[k]]
    long_range <- fit$coefficients$long_range[t1[k], t2[k]]
    gamma <- fit$coefficients$gamma[t1[k], t2[k]]

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
  rev <- lapply(seq_len(length(t1)), eval_k)
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

compute_xseq <- function(potentials, nsteps = 1e3) {
  max_distance <- 1.5 * max(max(sapply(potentials$short_range, function(y) Reduce("max", y))), Reduce("max", potentials$long_range))
  seq(from = 0, to = max_distance, length.out = nsteps)
}

format_potentials <- function(pot,
                              epsilon = 1e-4) { # epsilon controls the change in potential that is worth considering when classifying curves
  str <- paste0("Interaction potentials corresponding to a provided fit.\n\n",
                "The object contains three slots:\n- obj$short is a list containing functions of x, ",
                "each entry corresponding to one of the short-range potentials.\n",
                "- obj$medium contains a function of x representing the medium-range potential.\n",
                "- obj$overall contains a function of x representing the sum of all potentials.\n")

  is_one_pair <- length(pot$type1) == 1
  if(!is_one_pair) {
    xseq <- compute_xseq(pot)

    prs <- paste0(pot$type1, pot$type2)

    possible_behaviours <- c("Always decreasing",
                             "Always increasing",
                             "Decreasing, then increasing",
                             "Increasing, then decreasing")
    behaviour_pair <- setNames(sapply(unique(prs), function(pr) {
      val <- pot$overall[[which(unique(prs) == pr)]](xseq)
      val <- val[!is.na(val)]
      first_change <- val[-1][abs(val[-1] - val[1]) > epsilon][1]
      ifelse(all(sort(val, decreasing = TRUE) == val),
             possible_behaviours[1],
             ifelse(all(sort(val, decreasing = FALSE) == val),
                    possible_behaviours[2],
                    ifelse(first_change < val[1],
                           possible_behaviours[3],
                           possible_behaviours[4])))
    }), nm = unique(prs))

    behaviour_pair <- setNames(sapply(possible_behaviours, function(b) sum(behaviour_pair == b)), nm = possible_behaviours)
    str <- paste0(str, "\nSummary of curve behaviour: ", paste0(names(behaviour_pair), " = ", behaviour_pair,  collapse = "; "), ".\n")
  }
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
plot.potentials <- function(x, base_size = 11, max_npairs = 100, compute_average = TRUE, compute_median = TRUE, ...) {
  xseq <- compute_xseq(x)

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

  npairs <- length(unique(df$pair))

  if(npairs > 1 & (compute_average | compute_median)) {
    pairs <- unique(df$pair)

    avg <- lapply(pairs, function(pr) {
      df[df$pair == pr, c("x", "Overall")]
    })

    avg_df <- data.frame(x = unique(Reduce(c, lapply(avg, function(a) a$x))))

    if(compute_median) {
      avg_df$median <- sapply(avg_df$x, function(x) median(sapply(avg, function(a) a$Overall[a$x == x])))
    }

    if(compute_average) {
      avg_df$average <- sapply(avg_df$x, function(x) mean(sapply(avg, function(a) a$Overall[a$x == x])))
    }
  }

  g <- ggplot(data = df) +
    geom_line(aes(x = .data$x, y = 0), colour = "black") +
    xlab(NULL) +
    ylab(NULL) +
    theme_minimal(base_size = base_size) +
    theme(legend.title = element_blank())

  # Hide legend and make curves more transparent if too many pairs
  if(npairs > max_npairs) {
    alpha_main <- 0.3
    alpha_secondary <- 0.1
  } else {
    alpha_main <- 0.6
    alpha_secondary <- 0.5
  }

  if(npairs == 1) {
    g <- g + geom_line(aes(x = .data$x, y = .data$Overall, colour = "Overall", linetype = "Overall"), size = 2)
  } else {
    g <- g + geom_line(aes(x = .data$x, y = .data$Overall, linetype = "Overall", colour = .data$pair), size = 1.5, alpha = alpha_main)
  }


  nshorts <- sum(startsWith(colnames(df), "short"))
  for(i in seq_len(nshorts)) {
    if(npairs == 1) {
      g <- g + geom_line(aes(x = .data$x, y = .data[[paste0("short", i)]], linetype = .data[[paste0("name", i)]],
                             colour = .data[[paste0("name", i)]]), size = 1.5, alpha = 0.7)
    } else {
      g <- g + geom_line(aes(x = .data$x, y = .data[[paste0("short", i)]],
                             linetype = "Short/Medium", colour = .data$pair), size = 1, alpha = alpha_secondary)
    }
  }

  if(any(df$Medium != 0, na.rm = TRUE)) {
    if(npairs == 1) {
      g <- g + geom_line(aes(x = .data$x, y = .data$Medium, colour = "Medium", linetype = "Medium"), size = 1.5, alpha = 0.7)
    } else {
      g <- g + geom_line(aes(x = .data$x, y = .data$Medium, linetype = "Short/Medium", colour = .data$pair), size = 1, alpha = alpha_secondary)
    }
  }

  if(npairs == 1) {
    lab <- "Overall"
    for(i in seq_len(nshorts)) {
      lab[length(lab) + 1] <- df[1, paste0("name", i)]
    }
    if(any(df$Medium != 0, na.rm = TRUE)) {
      lab[length(lab) + 1] <- "Medium"
    }
    val <- "solid"
    if(length(lab) > 1) {
      val <- c(val, rep(c("dotted", "longdash"), length(lab) - 1))
    }
    val <- val[order(lab)]
    lab <- sort(lab)
    g <- g +
      scale_colour_viridis_d(name = "", labels = lab, begin = 0.1, end = 0.9, option = "turbo") +
      scale_linetype_manual(name = "", values = val, labels = lab) +
      guides(colour = guide_legend(keywidth = 5, override.aes = list(linewidth = 1.5)))
  } else {
    if(npairs > max_npairs) {
      g <- g + scale_colour_viridis_d(guide = "none")
    } else {
      g <- g + scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "turbo",
                                      guide = guide_legend(order = 2, nrow = 10))
    }
    g <- g +
      scale_linetype_manual(values = c("solid", "dotted", "longdash")) +
      guides(linetype = guide_legend(order = 3,
                                     keywidth = 3,
                                     override.aes = list(linewidth = 1))) +
      new_scale_colour()
    if(compute_median) {
      g <- g + geom_line(data = avg_df, aes(x = .data$x, y = .data$median, colour = "Median"), alpha = 0.9, size = 2.5)
    }
    if(compute_average) {
      g <- g + geom_line(data = avg_df, aes(x = .data$x, y = .data$average, colour = "Average"), alpha = 0.9, size = 2.5)
    }
    g <- g + scale_colour_manual(values = c(Average = "black", Median = "red"), guide = guide_legend(order = 1))
  }

  g
}
