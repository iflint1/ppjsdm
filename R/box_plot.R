#' Construct a list of fits from a series of fit objects.
#'
#' @param ... Any number of fit objects obtained by a call to `ppjsdm::gibbsm`.
#' @param list Some more fits provided as a list.
#' @export
#' @keywords internal
#' @md
create_list_of_fits <- function(...,
                                list) {
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

  # If user did not give names to the fits, add some default names
  default_names <- paste0("Fit ", seq_len(length(fits)))
  if(is.null(names(fits))) {
    names(fits) <- default_names
  }
  names(fits) <- ifelse(names(fits) == "", default_names, names(fits))

  fits
}

#' Starting from a list of fits, provide a function to access a coefficient given by the user as a string.
#'
#' @param fits A list of fit objected obtained by a call to `ppjsdm::gibbsm`.
#' @param coefficient A string representing a coefficient of interest.
#' Choice of `alpha1`, `alpha2`, ..., to show one of the short-range interaction coefficients;
#' `gamma` to show the medium-range interaction coefficient;
#' `beta1`, `beta2`, ... or actual name of the covariate for one of the regression coefficients.
#' `alpha` shows a facet plot of all short-range interaction coefficients, while `beta` does the same for the regression coefficients.
#' @export
#' @keywords internal
#' @md
access_coefficient <- function(fits,
                               coefficient) {
  if("gamma" == coefficient[1] & length(coefficient) == 1) { # The requested coefficient is gamma
    identification <- "gamma"
    access <- function(obj) base::list(obj[["gamma"]])
    index <- NULL
  } else {
    is_alpha <- gsub("^alpha([0-9]*)", "\\1", coefficient)
    if("" == is_alpha[1] & length(is_alpha) == 1) { # it starts with alpha, but has no integer afterwards
      index <- seq_len(max(sapply(fits, function(fit) {
        length(fit$coefficients$alpha)
      }))) # all the alphas are of interest
      access <- function(obj) obj[["alpha"]]
      identification <- "alpha"
    } else if(length(coefficient) == 1 & is_alpha[1] != coefficient[1]) { # it starts with alpha, with an integer afterwards
      identification <- "alpha"
      access <- function(obj) obj[["alpha"]][as.numeric(is_alpha)]
      index <- is_alpha
    } else {
      is_beta <- gsub("^beta([0-9]*)", "\\1", coefficient)
      identification <- "beta"
      list_covariates <- lapply(fits, function(fit) { # Contains all the covariates
        colnames(as.matrix(fit$coefficients[["beta"]]))
      })
      list_covariates <- as.character(unique(Reduce(c, list_covariates)))
      if("" == is_beta[1] & length(is_beta) == 1) { # it starts with beta, but has no integer afterwards
        index <- list_covariates
      } else if(!all(is.na(suppressWarnings(as.numeric(is_beta))))) { # is_beta is an integer, telling us which beta to access
        index <- as.numeric(is_beta)
      } else if(all(coefficient %in% list_covariates)) {
        index <- Reduce(c, coefficient)
      } else {
        stop("Unrecognised format for the coefficient.")
      }

      # If it starts with beta, then this works. If not, assume that coefficient is a covariate and access in this way too.
      access <- function(obj) {
        beta <- as.matrix(obj[["beta"]])
        if(is.character(index)) {
          i <- intersect(index, colnames(beta))
          z <- as.matrix(beta[, i])
          colnames(z) <- i
          base::list(z)
        } else {
          base::list(beta)
        }
      }
    }
  }

  list(access = access,
       identification = identification,
       index = index)
}

#' Convert a vector of names to their corresponding entry in a named list.
#'
#' @param x A vector of names.
#' @param named_list A named list, where names correspond to EITHER elements in x OR x converted according to one of the extra supplied objects
#' , and the value in `named_list` represent what we are converting to.
#' @param ... We assume that the names in `named_list` might correspond to a conversion of x according to named lists in these objects.
#' @export
#' @keywords internal
#' @md
convert_names <- function(x,
                          named_list,
                          ...) {
  possible_names <- list(...)
  if(!is.list(named_list)) {
    stop("named_list should be a named list, with names corresponding to the abbreviated names.")
  }
  sapply(as.character(x), function(ty) {
    if(length(named_list[[ty]]) > 0) {
      return(named_list[[ty]])
    }
    for(possible_names in possible_names) {
      if(!is.list(possible_names)) {
        stop("The extra supplied objects should all be (named) lists.")
      }
      if(length(possible_names[[ty]]) > 0) {
        if(length(named_list[[possible_names[[ty]]]]) > 0) {
          return(named_list[[possible_names[[ty]]]])
        }
      }
    }
    ty
  })
}

#' Construct a `data.frame` containing estimated coefficients and their confidence intervals from some fit objects.
#'
#' This function is used to prepare data for use in the `box_plot` and `diagram_plot` functions.
#'
#' @param fits A list of fit objected obtained by a call to `ppjsdm::gibbsm`.
#' @param coefficient A string representing the coefficient to plot.
#' Choice of `alpha1`, `alpha2`, ..., to show one of the short-range interaction coefficients;
#' `gamma` to show the medium-range interaction coefficient;
#' `beta1`, `beta2`, ... or actual name of the covariate for one of the regression coefficients.
#' `alpha` shows a facet plot of all short-range interaction coefficients, while `beta` does the same for the regression coefficients.
#' @param summ Optional list of summaries corresponding to the fits; if not provided they are obtained by calling `ppjsdm::summary.gibbsm`
#' @param only_statistically_significant Only show statistically significant coefficients?
#' @param which If plotting interaction coefficients, which ones do we want to plot?
#' @param full_names Optional list of full names of types, if for example abbreviations were used when running the fit.
#' @param compute_confidence_intervals Compute the confidence intervals (which is slower) or just show the point estimates?
#' @param classes If this parameter is supplied, then colours are used to distinguish classes instead of fits.
#' Only works when given a single fit.
#' Should be a named vector/list, with the names corresponding to types, and the value equal to the class name.
#' @param involving Optional vector/list of types. Only coefficients involving these types will be plotted.
#' @param how If the `involving` argument is supplied, should it involve *only* those types, or at least *one* of those types (which is relevant if inter-type interactions are involved).
#' @export
#' @keywords internal
#' @md
make_summary_df <- function(fits,
                            coefficient,
                            summ,
                            only_statistically_significant,
                            which,
                            full_names,
                            compute_confidence_intervals,
                            classes,
                            involving,
                            how) {
  # Take care of the classes argument
  if(!is.null(classes)) {
    classes <- as.list(classes)
    if(length(fits) > 1) {
      stop("If classes is supplied, then there should be a single fit. Otherwise, colours are used to distinguish fits, not classes.")
    }
  }

  # Take care of the full_names argument
  if(!is.null(full_names)) {
    full_names <- as.list(full_names)
    if(length(unique(full_names)) != length(full_names)) {
      stop("full_names does not contain unique names: this will cause some issues later on, please supply distinct names.")
    }
  }

  read_coefficient <- access_coefficient(fits = fits, coefficient = coefficient)
  access <- read_coefficient$access
  identification <- read_coefficient$identification
  index <- read_coefficient$index

  # If coefficient is not an interaction parameter, need to do extra work
  if(identification == "beta" & which != "all") {
    warning("Plotting regression coefficient but parameter \"which\" was set to something other than \"all\". Assuming this is a typo and setting \"which\" to \"all\".")
    which <- "all"
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

  # Since our estimates are a list of fits, convert df to a list of dataframes, each one corresponding to one of the fits.
  dfs <- lapply(seq_len(length(estimates)), function(k) {
    dfs_by_potential <- lapply(seq_len(length(estimates[[k]])), function(n) {
      # Get the types of the current fit
      types <- rownames(estimates[[k]][[n]])

      # Types involving the focal ones
      types_subset <- if(!is.null(involving)) {
        involving <- Reduce(c, as.list(involving))
        intersect(types, involving)
      } else {
        types
      }

      # Should other types be considered?
      other_types_subset <- if(how == "only") {
        types_subset
      } else {
        types
      }

      if(identification != "beta") { # Create dataframe of two columns with all possible pairs of types
        d <- as.data.frame(expand.grid(from = other_types_subset, to = types_subset, stringsAsFactors = FALSE))
        d <- d[!duplicated(t(apply(d, 1, sort))), ]

        if(which == "within") { # in this case, remove columns that are the same, so shows only between interactions
          d <- d[d$to == d$from, ]
        } else if(which == "between" | which == "average_between") {
          d <- d[d$to != d$from, ]
        }
      } else if(is.numeric(index)) {
        d <- as.data.frame(Reduce(rbind, lapply(index, function(i) data.frame(from = types_subset, to = i))))
      } else {
        d <- as.data.frame(Reduce(rbind, lapply(intersect(index, colnames(fits[[k]]$coefficients$beta)), function(i) data.frame(from = types_subset, to = i))))
      }

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

        d$lo_numerical <- sapply(seq_len(nrow(d)), function(i) { # Lower-endpoint of the numerical CI (relating to the numerical error due to dummy points)
          tryCatch(access(summ[[k]]$lo_numerical)[[n]][d$from[i], d$to[i]], error = function(err) NA)
        })

        d$hi_numerical <- sapply(seq_len(nrow(d)), function(i) { # Upper-endpoint of the numerical CI (relating to the numerical error due to dummy points)
          tryCatch(access(summ[[k]]$hi_numerical)[[n]][d$from[i], d$to[i]], error = function(err) NA)
        })
      }

      # Convert names to full names
      if(!is.null(full_names)) {
        # We do not yet convert the names, because to decide class we want to look at original names first
        d$new_from <- convert_names(d$from, full_names)
        if(identification != "beta") {
          d$new_to <- convert_names(d$to, full_names)
        }
      }

      # If classes was supplied, fill in class
      if(!is.null(classes)) {
        if(!is.null(full_names)) {
          d$class_from <- convert_names(d$from, classes, full_names)
          d$class_to <- convert_names(d$to, classes, full_names)
        } else {
          d$class_from <- convert_names(d$from, classes)
          d$class_to <- convert_names(d$to, classes)
        }
      }

      # At this point, we can discard the temporary column names chosen above
      if("new_from" %in% colnames(d)) {
        d$from <- d$new_from
        d$new_from <- NULL
      }
      if("new_to" %in% colnames(d)) {
        d$to <- d$new_to
        d$new_to <- NULL
      }

      # Set name of the fit
      if(nrow(d) > 0) {
        d$Fit <- names(fits)[k]
        d$Potential <- paste0("Potential ", index[n])
      }

      # Compute average of the coefficient for each type
      if(which == "average_between") {
        d <- as.data.frame(Reduce(rbind, lapply(union(unique(d$from), unique(d$to)), function(ty) {
          g <- d[d$from == ty | d$to == ty, ][1, ]
          g$from <- ty
          g$to <- ty
          g$E <- mean(d$E[d$from == ty | d$to == ty], na.rm = TRUE)
          g
        })))
      }

      if(only_statistically_significant) { # in this case, remove non-stat significant
        d <- as.data.frame(d[d$lo > 0 | d$hi < 0, ]) # If non-statistically significant, remove row,. i.e. if low is > 0 or if high < 0, the row is kept
      }

      d
    })

    Reduce(rbind, dfs_by_potential)
  })

  # Flatten dfs
  df <- as.data.frame(Reduce(rbind, dfs))

  # Make into factor and order them correctly
  if(!is.null(df$Fit)) {
    df$Fit <- factor(df$Fit, levels = names(fits))
  }
  if(!is.null(df$Potential)) {
    df$Potential <- factor(df$Potential)
  }
  if(!is.null(df$class_from)) {
    df$class_from <- factor(df$class_from)
  }
  if(!is.null(df$class_to)) {
    df$class_to <- factor(df$class_to)
  }

  colnames(df)[colnames(df) == "E"] <- identification

  df
}

#' Plot the coefficients of fit objects.
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
#' @param summ Optional list of summaries corresponding to the fits; if not provided they are obtained by calling `ppjsdm::summary.gibbsm`
#' @param only_statistically_significant Only show statistically significant coefficients?
#' @param which If plotting interaction coefficients, which ones do we want to plot?
#' @param full_names Optional list of full names of types, if for example abbreviations were used when running the fit.
#' @param compute_confidence_intervals Compute the confidence intervals (which is slower) or just show the point estimates?
#' @param classes If this parameter is supplied, then colours are used to distinguish classes instead of fits.
#' Only works when given a single fit.
#' Should be a named vector/list, with the names corresponding to types, and the value equal to the class name.
#' @param involving Optional vector/list of types. Only coefficients involving these types will be plotted.
#' @param how If the `involving` argument is supplied, should it involve *only* those types, or at least *one* of those types (which is relevant if inter-type interactions are involved).
#' @param title Plot title.
#' @param colours Optional vector of colours to represent the different fits/classes.
#' @param highlight_zero Highlight the zero value with a red line?
#' @param text_size Text size.
#' @param base_size Base size.
#' @param xmin Optional plot minimum x.
#' @param xmax Optional plot maximum x.
#' @importFrom ggplot2 aes element_text facet_grid geom_errorbar geom_point geom_vline ggplot ggtitle guide_legend guides position_dodge scale_color_manual theme theme_bw xlab ylab xlim
#' @importFrom rlang .data
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
#' ppjsdm::box_plot(fit)
#'
#' ppjsdm::box_plot(fit, coefficient = "x")
#'
#' @export
#' @md
box_plot <- function(...,
                     list,
                     coefficient = "alpha1",
                     summ,
                     only_statistically_significant = FALSE,
                     which = c("all", "within", "between", "average_between"),
                     full_names = NULL,
                     compute_confidence_intervals = TRUE,
                     classes = NULL,
                     involving = NULL,
                     how = c("only", "one"),
                     title,
                     colours = c("black", "#5BBCD6", "#F2AD00", "#00A08A", "#FF0000"),
                     highlight_zero = TRUE,
                     text_size = 16,
                     base_size = 20,
                     xmin,
                     xmax) {
  # Interpret the how argument
  how <- match.arg(how)

  # TODO: Avoid NULL default arguments for args forwarded to internal function, just set those to NULL in there
  # TODO: Is there a clean way to copy the documentation between functions when they have the same arguments? See box_plot and chord_diagram, all the gibbsm rgibbs functions, etc.
  # Interpret the which argument
  which <- match.arg(which)

  fits <- create_list_of_fits(..., list = list)

  df <- make_summary_df(fits = fits,
                        coefficient = coefficient,
                        summ = summ,
                        only_statistically_significant = only_statistically_significant,
                        which = which,
                        full_names = full_names,
                        compute_confidence_intervals = compute_confidence_intervals,
                        classes = classes,
                        involving = involving,
                        how = how)

  # Name of the coefficient we are plotting
  identification <- c("gamma", "alpha", "beta")
  identification <- identification[identification %in% colnames(df)]

  # Remove NAs
  df <- df[!is.na(df[, identification]), ]

  if(length(identification) != 1) {
    stop("Could not identify the intended coefficient by analysing the colnames of the dataframe.")
  }

  # Make default title
  if(missing(title)) {
    if(identification == "gamma") {
      title <- "Medium-range interaction coefficients"
    } else if(identification == "alpha") {
      if(length(unique(df$Potential)) > 1) {
        title <- "Short-range interaction coefficients"
      } else {
        title <- paste0("Short-range interaction coefficients for ", df$Potential)
      }
    } else if(identification == "beta") {
      if(length(unique(df$to)) > 1) {
        title <- "Beta coefficients"
      } else if(is.numeric(df$to)) {
        title <- paste0("Beta coefficients for covariate number ", df$to)
      } else {
        title <- paste0("Beta coefficients for ", df$to)
      }
    } else {
      stop("Cannot come up with a good title; did not recognise the coefficient being plotted.")
    }
  }

  # Define and order y-axis labels
  nc <- nrow(df)
  ylabels <- factor(seq_len(nc))
  if(identification != "beta" && which != "within" && which != "average_between") {
    levels(ylabels) <- ifelse(df$from < df$to, paste0(df$from, enc2utf8(" \u2194 "), df$to), paste0(df$to, enc2utf8(" \u2194 "), df$from)) # Order lhs and rhs of interaction alphabetically
  } else {
    levels(ylabels) <- df$from
  }
  df$ylabels <- ylabels

  # Do we want a facet plot? How about ordering of y labels?
  if(identification == "beta" & length(unique(df$to)) > 1) { # In this case, alphabetical order + facet plot
    do_split <- TRUE
    df$split <- df$to
    df$ylabels <- factor(df$ylabels, levels = sort(levels(df$ylabels), decreasing = TRUE))
  } else if(identification == "alpha" & length(unique(df$Potential)) > 1) { # Same as above, but facet plot on potential name
    do_split <- TRUE
    df$split <- df$Potential
    df$ylabels <- factor(df$ylabels, levels = sort(levels(df$ylabels), decreasing = TRUE))
  } else {
    # Compute average (over the fits) estimate for each coefficient
    average_over_fits <- sapply(levels(df$ylabels), function(ylab) {
      mean(df[df$ylabels == ylab, identification], na.rm = TRUE)
    })

    df$ylabels <- factor(df$ylabels, levels = levels(df$ylabels)[order(average_over_fits)])
    do_split <- FALSE
  }

  # What do colours represent?
  if(is.null(df$class_from)) {
    colour_string <- "Fit"
    ncolours <- nlevels(df$Fit)
  } else {
    colour_string <- "Class"
    df$Class <- df$class_from
    ncolours <- nlevels(df$class_from)
  }

  # Set colours
  colours <- rep(colours, length.out = ncolours)

  g <- ggplot(df, aes(y = .data$ylabels, x = .data[[identification]])) # Visualisation

  if(highlight_zero) { # If this option is chosen, draw a red line before doing anything else
    g <- g + geom_vline(xintercept = 0, color = "red", linewidth = 1.5) # Vertical line at 0
  }

  g <- g +
    geom_point(aes(colour = .data[[colour_string]]), size = 3, position = position_dodge(width = 0.2)) + # Size of point estimates
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

  if(all(c("lo_numerical", "hi_numerical", "lo", "hi") %in% colnames(df))) {
    g <- g + geom_errorbar(aes(colour = .data[[colour_string]],
                               xmin = .data$lo_numerical,
                               xmax = .data$hi_numerical), width = 0, linewidth = 2, position = position_dodge(width = 0.2)) +
      geom_errorbar(aes(colour = .data[[colour_string]],
                        xmin = .data$lo,
                        xmax = .data$hi), width = 0.2, linewidth = 0.5, position = position_dodge(width = 0.2))
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
