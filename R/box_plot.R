#' Plot the coefficients of a fit object.
#'
#' The coefficients and corresponding confidence intervals are shown in a box-like plot.
#' The inner thick parts of the error bars represent numerical uncertainty, due to the number and distribution of dummy points.
#' The outer thin part of the error bars represent the total statistical + numerical uncertainty.
#' Formally, these are theoretical asymptotic 95% confidence intervals.
#' As a convenience, this function can alternatively be applied by calling `plot` on a fit object.
#'
#' @param ... Any number of fit objects obtained by a call to `ppjsdm::gibbsm`.
#' @param list Some more fits provided as a list.
#' @param coefficient A string representing the coefficient to plot.
#' Choice of `alpha1`, `alpha2`, ..., to show one of the short-range interaction coefficients;
#' `gamma` to show the medium-range interaction coefficient;
#' `beta1`, `beta2`, ... or actual name of the covariate for one of the regression coefficients.
#' `alpha` shows a facet plot of all short-range interaction coefficients, while `beta` does the same for the regression coefficients.
#' @param title Plot title.
#' @param summ Optional list of summaries corresponding to the fits; if not provided they are obtained by calling `ppjsdm::summary.gibbsm`
#' @param only_statistically_significant Only show statistically significant coefficients?
#' @param which If plotting interaction coefficients, which ones do we want to plot?
#' @param highlight_zero Highlight the zero value with a red line?
#' @param text_size Text size.
#' @param base_size Base size.
#' @param full_names Optional list of full names of types, if for example abbreviations were used when running the fit.
#' @param compute_confidence_intervals Compute the confidence intervals (which is slower) or just show the point estimates?
#' @param colours Optional vector of colours to represent the different fits.
#' @param classes If this parameter is supplied, then colours are used to distinguish classes instead of fits.
#' Only works when given a single fit.
#' Should be a named vector/list, with the names corresponding to types, and the value equal to the class name.
#' @param involving Optional vector/list of types. Only coefficients involving these types will be plotted.
#' @param xmin Optional plot minimum x.
#' @param xmax Optional plot maximum x.
#' @importFrom ggplot2 aes aes_string element_text facet_grid geom_errorbar geom_point geom_vline ggplot ggtitle guide_legend guides position_dodge scale_color_manual theme theme_bw xlab ylab xlim
#' @examples
#' set.seed(1)
#'
#' # Construct a configuration
#'
#' configuration <- ppjsdm::rppp(lambda = c(A = 500, B = 500))
#'
#' # Fit the model
#'
#' fit <- ppjsdm::gibbsm(configuration, covariates = list(x = function(x, y) x))
#'
#' # Plot the coefficients
#'
#' box_plot(fit)
#'
#' box_plot(fit, coefficient = "x")
#'
#' @export
#' @md
box_plot <- function(...,
                     list,
                     coefficient = "alpha1",
                     title,
                     summ,
                     only_statistically_significant = FALSE,
                     which = c("all", "within", "between", "average_between"),
                     highlight_zero = TRUE,
                     text_size = 16,
                     base_size = 20,
                     full_names = NULL,
                     compute_confidence_intervals = TRUE,
                     colours,
                     classes = NULL,
                     involving,
                     xmin,
                     xmax) {
  # Interpret the which argument
  which <- match.arg(which)

  # Allow for either sequence of fits or list of fits, convert both to list
  fits <- if(missing(list)) {
    base::list(...)
  } else {
    c(base::list(...), list)
  }

  # Subset to objects which have the right class, discarding others
  which_gibbsm <- sapply(fits, function(fit) inherits(fit, "gibbsm"))
  if(!all(which_gibbsm)) {
    warning("Provided some arguments that were not fit objects, these were discarded.")
  }
  if(length(which_gibbsm) == 0) {
    stop("No valid fits were provided.")
  }
  fits <- fits[which_gibbsm]

  # Take care of the classes argument
  if(!is.null(classes)) {
    classes <- as.list(classes)
    if(length(fits) > 1) {
      stop("If classes is supplied, then there should be a single fit. Otherwise, colours are used to distinguish fits, not classes.")
    }
  }

  # If user did not give names to the fits, add some default names
  default_names <- paste0("Fit ", seq_len(length(fits)))
  if(is.null(names(fits))) {
    names(fits) <- default_names
  }
  names(fits) <- ifelse(names(fits) == "", default_names, names(fits))

  # Below, we translate the coefficient parameter into a function "access" that allows us to access it in the fit object
  # In addition, is_interaction tells us whether or not the coefficient is an interaction coefficient (alpha/gamma) or not (beta)
  is_interaction <- TRUE
  is_alpha <- NULL
  is_beta <- NULL
  if("gamma" == coefficient[1] & length(coefficient) == 1) { # The requested coefficient is gamma
    if(missing(title)) {
      title <- "Medium-range interaction coefficients"
    }
    access <- function(obj) base::list(obj[["gamma"]])
  } else {
    is_alpha <- gsub("^alpha([0-9]*)", "\\1", coefficient)
    if("" == is_alpha[1] & length(is_alpha) == 1) { # it starts with alpha, but has no integer afterwards
      is_alpha <- seq_len(max(sapply(fits, function(fit) {
        length(fit$coefficients$alpha)
      }))) # show a facet plot of all potentials
      if(missing(title)) {
        title <- "Short-range interaction coefficients"
      }
    }

    if(length(coefficient) == 1 & is_alpha[1] != coefficient[1]) { # it starts with alpha
      if(missing(title)) {
        title <- paste0("Short-range interaction coefficients for potential ", is_alpha[1])
      }
      access <- function(obj) obj[["alpha"]][as.numeric(is_alpha)]
    } else {
      is_interaction <- FALSE
      if(which != "") {
        warning("Plotting regression coefficient but parameter \"which\" was set to something other than \"all\". Assuming this is a typo and setting \"which\" to \"all\".")
        which <- "all"
      }
      is_beta <- gsub("^beta([0-9]*)", "\\1", coefficient)
      if("" == is_beta[1] & length(is_beta) == 1) { # it starts with beta, but has no integer afterwards
        # In this case, we want to get all covariates involved, and use all
        list_covariates <- lapply(fits, function(fit) {
          colnames(as.matrix(fit$coefficients[["beta"]]))
        })
        is_beta <- as.character(unique(Reduce(c, list_covariates)))
        if(missing(title)) {
          title <- "Beta coefficients"
        }
      } else if(!all(is.na(suppressWarnings(as.numeric(is_beta))))) { # is_beta is an integer, telling us which beta to access
        is_beta <- as.numeric(is_beta)
        if(missing(title)) {
          title <- paste0("Beta coefficient for covariate number ", is_beta)
        }
      } else if(missing(title)) {
        title <- paste0("Beta coefficient for ", is_beta)
      }

      # If it starts with beta, then this works. If not, assume that coefficient is a covariate and access in this way too.
      access <- function(obj) {
        beta <- as.matrix(obj[["beta"]])
        if(is.character(is_beta)) {
          i <- intersect(is_beta, colnames(beta))
          z <- as.matrix(beta[, i])
          colnames(z) <- i
          base::list(z)
        } else {
          base::list(beta)
        }
      }
    }
  }

  # Do not compute CIs if average_between
  compute_confidence_intervals <- compute_confidence_intervals & which != "average_between"

  if(compute_confidence_intervals) {
    if(!missing(summ)) { # Make sure summaries and fits are compatible
      if(!is(summ, "list")) {
        summ <- base::list(summ)
      }
      stopifnot(length(fits) == length(summ))
    } else { # Construct the summaries
      summ <- lapply(fits, function(f) summary(f))
    }
  }

  # Extract list of coefficients, each one corresponding to one of the fits
  estimates <- lapply(fits, function(f) {
    tryCatch(access(f$coefficients), error = function(err) NA)
  })

  # The fits could perhaps relate to fits which do not involve the same types.
  # Construct a set of all types that appear in some of the fits
  types <- lapply(estimates, function(e) unique(Reduce(c, lapply(e, function(f) rownames(f)))))
  types <- unique(Reduce(c, types))

  # Types involving the focal ones
  types_subset <- if(!missing(involving)) {
    intersect(types, involving)
  } else {
    types
  }

  if(is_interaction) { # Create dataframe of two columns with all possible pairs of types
    df <- as.data.frame(expand.grid(from = types, to = types_subset, stringsAsFactors = FALSE))
    df <- df[!duplicated(t(apply(df, 1, sort))), ]

    if(which == "within") { # in this case, remove columns that are the same, so shows only between interactions
      df <- df[df$to == df$from, ]
    } else if(which == "between" | which == "average_between") {
      df <- df[df$to != df$from, ]
    }
  } else {
    df <- Reduce(rbind, lapply(is_beta, function(i) data.frame(from = types_subset, to = i)))
  }

  # At this point, df is a single dataframe containing all relevant pairs of types.
  # Since our estimates are a list of fits, convert df to a list of dataframe, each one corresponding to one of the fits.
  dfs <- lapply(seq_len(length(estimates)), function(k) {
    dfs_by_potential <- lapply(seq_len(length(estimates[[k]])), function(n) {
      d <- df

      d$E <- sapply(seq_len(nrow(d)), function(i) { #these and below parts specify the values of point and bars
        tryCatch(estimates[[k]][[n]][d$from[i], d$to[i]], error = function(err) NA)
      })

      if(compute_confidence_intervals) {
        d$lo <- sapply(seq_len(nrow(d)), function(i) { # Get the lower-endpoint of the CIs, apply function to sequence of numbers 1 to nrow.df
          tryCatch(access(summ[[k]]$lo)[[n]][d$from[i], d$to[i]], error = function(err) NA) # retrieve numbers from summ and apply to the rows
        })

        d$hi <- sapply(seq_len(nrow(d)), function(i) { # Get the upper-endpoint of the CIs
          tryCatch(access(summ[[k]]$hi)[[n]][d$from[i], d$to[i]], error = function(err) NA)
        })

        if(only_statistically_significant) { # in this case, remove non-stat significant
          d <- d[d$lo > 0 | d$hi < 0, ] # If non-statistically significant, remove row,. i.e. if low is > 0 or if high < 0, the row is kept
        }

        d$lo_numerical <- sapply(seq_len(nrow(d)), function(i) { # Lower-endpoint of the numerical CI (relating to the numerical error due to dummy points)
          tryCatch(access(summ[[k]]$lo_numerical)[[n]][d$from[i], d$to[i]], error = function(err) NA)
        })

        d$hi_numerical <- sapply(seq_len(nrow(d)), function(i) { # Upper-endpoint of the numerical CI (relating to the numerical error due to dummy points)
          tryCatch(access(summ[[k]]$hi_numerical)[[n]][d$from[i], d$to[i]], error = function(err) NA)
        })
      }

      # If classes was supplied, fill in class
      if(!is.null(classes)) {
        d$Class <- sapply(as.character(d$from), function(sp)  classes[[sp]])
      }

      if(!is.null(full_names)) { # specifying names
        full_names <- as.list(full_names)

        d$from <- sapply(as.character(d$from), function(sp) if(length(full_names[[sp]]) == 0) sp else full_names[[sp]])
        if(is_interaction) {
          d$to <- sapply(as.character(d$to), function(sp) if(length(full_names[[sp]]) == 0) sp else full_names[[sp]])
        }
      }

      # Set name of the fit
      d$Fit <- names(fits)[k]
      d$Potential <- paste0("Potential ", n)

      # Compute average of the coefficient for each type
      if(which == "average_between") {
        d <- Reduce(rbind, lapply(union(unique(d$from), unique(d$to)), function(ty) {
          g <- d[d$from == ty | d$to == ty, ][1, ]
          g$from <- ty
          g$to <- ty
          g$E <- mean(d$E[d$from == ty | d$to == ty], na.rm = TRUE)
          g
        }))
      }

      d
    })

    Reduce(rbind, dfs_by_potential)
  })

  # Flatten dfs
  df <- Reduce(rbind, dfs)

  # Compute average (over the fits) estimate for each coefficient
  df$average_estimates <- sapply(seq_len(nrow(df)), function(i) {
    mean(df$E[df$from == df$from[i] & df$to == df$to[i]], na.rm = TRUE)
  })

  # Define and order y-axis labels
  nc <- nrow(df)
  x <- factor(1:nc)
  if(is_interaction && which != "within" && which != "average_between") {
    levels(x) <- ifelse(df$from < df$to, paste0(df$from, enc2utf8(" \u2194 "), df$to), paste0(df$to, enc2utf8(" \u2194 "), df$from)) # Order lhs and rhs of interaction alphabetically
  } else {
    levels(x) <- df$from
  }

  if(!is_interaction  & length(is_beta) > 1) { # In this case, alphabetical order + facet plot
    do_split <- TRUE
    df$split <- df$to
    df$x <- factor(x, levels = sort(levels(x), decreasing = TRUE))
  } else if(is_interaction & length(is_alpha) > 1) { # Same as above, but facet plot on potential name
    do_split <- TRUE
    df$split <- df$Potential
    df$x <- factor(x, levels = sort(levels(x), decreasing = TRUE))
  } else {
    df$x <- factor(x, levels = levels(x)[order(df$average_estimates)])
    do_split <- FALSE
  }

  # Order Fit legend according to fits object
  df$Fit <- factor(df$Fit, levels = names(fits))

  # Remove NAs
  df <- df[rowSums(is.na(df)) == 0, ]

  # What do colours represent?
  if(is.null(classes)) {
    colour_string <- "Fit"
    ncolours <- nlevels(df$Fit)
  } else {
    df$Class <- factor(df$Class)
    colour_string <- "Class"
    ncolours <- nlevels(df$Class)
  }

  # Set colours
  if(missing(colours)) {
    colours <- rep(c("black", "#5BBCD6", "#F2AD00", "#00A08A", "#FF0000"), length.out = ncolours)
  } else {
    colours <- rep(colours, length.out = ncolours)
  }

  g <- ggplot(df, aes_string(y = "x", x = "E")) # Visualisation

  if(highlight_zero) { # If this option is chosen, draw a red line before doing anything else
    g <- g + geom_vline(xintercept = 0, color = "red", linewidth = 1.5) # Vertical line at 0
  }

  g <- g +
    geom_point(aes_string(colour = colour_string), size = 3, position = position_dodge(width = 0.2)) + # Size of point estimates
    scale_color_manual(values = colours) + # Fit colours
    xlab(NULL) + # Remove x labels
    ylab(NULL) + # Remove y labels
    ggtitle(title) + # Title
    theme_bw(base_size = base_size) + # Theme
    theme(plot.title = element_text(size = text_size, hjust = 0.5),
          axis.text.y = element_text(size = text_size)) +
    guides(colour = guide_legend(title = "",
                                 keywidth = 2.5,
                                 override.aes = aes(size = 4, linewidth = 1)))

  if(compute_confidence_intervals) {
    g <- g + geom_errorbar(aes_string(colour = colour_string,
                                      xmin = "lo_numerical",
                                      xmax = "hi_numerical"), width = 0, linewidth = 2, position = position_dodge(width = 0.2)) +
      geom_errorbar(aes_string(colour = colour_string,
                               xmin = "lo",
                               xmax = "hi"), width = 0.2, linewidth = 0.5, position = position_dodge(width = 0.2))
  }

  # If xmin and xmax are given use those instead of defaults
  if(!missing(xmin) && !missing(xmax)) {
    g <- g + xlim(xmin, xmax) # Fix the x-axis range
  }

  # If there is a single fit, do not show legend
  if(colour_string == "Fit" & nlevels(df$Fit) == 1 | colour_string == "Class" & nlevels(df$Class) == 1) {
    g <- g + guides(colour = "none")
  }
  if(do_split) {
    g <- g + facet_grid(split ~ .)
  }
  g
}

#' Plot a `ppjsdm::gibbsm` object.
#'
#' This calls `ppjsdm::box_plot`, see that function's documentation for details.
#' @method plot gibbsm
#' @param x A fit object obtained by calling `ppjsdm::gibbsm`.
#' @param ... Forwarded to `ppjsdm::box_plot`.
#' @export
#' @md
plot.gibbsm <- function(x, ...) {
  box_plot(x, ...)
}
