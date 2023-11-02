#' Plot the coefficients of a fit object.
#'
#' We use a box-like plot.
#' The inner thick parts of the error bars represent numerical uncertainty, due to the number and distribution of dummy points.
#' The outer thin part of the error bars represent the total statistical + numerical uncertainty.
#' Formally, these are theoretical asymptotic 95% confidence intervals.
#' As a convenience, this function can alternatively be applied by calling `plot` on a fit object.
#'
#' @param ... Any number of fit objects obtained by a call to `ppjsdm::gibbsm`.
#' @param list Alternatively, fits provided as a list object.
#' @param coefficient A string representing the coefficient to plot.
#' Choice of alpha1, alpha2, ..., for one of the short-range interaction coefficients;
#' gamma for the medium-range interaction coefficient;
#' beta1, beta2, ... or actual name of the env. covariate for one of the regression coefficients.
#' @param title Plot title.
#' @param summ Optional list of summaries corresponding to the fits; if not provided they are obtained by calling `ppjsdm::summary.gibbsm`
#' @param only_statistically_significant Only show statistically significant coefficients?
#' @param which If plotting interaction coefficients, which ones do we want to plot?
#' @param highlight_zero Highlight the zero value with a red line?
#' @param text_size Text size.
#' @param base_size Base size.
#' @param full_names Optional list of full names of types, if for example abbreviations were used when running the fit.
#' @param colours Optional vector of colours to represent the different fits.
#' @param involving Optional vector/list of types. Only coefficients involving these types will be plotted.
#' @param xmin Optional plot minimum x.
#' @param xmax Optional plot maximum x.
#' @importFrom ggplot2 aes aes_string element_text facet_grid geom_errorbar geom_point geom_vline ggplot ggtitle guide_legend guides position_dodge scale_color_manual theme theme_bw xlab ylab xlim
#' @examples
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
                     title = "",
                     summ,
                     only_statistically_significant = FALSE,
                     which = c("all", "within", "between"),
                     highlight_zero = TRUE,
                     text_size = 16,
                     base_size = 20,
                     full_names = NULL,
                     colours,
                     involving,
                     xmin,
                     xmax) {
  # TODO: Add default title.
  # TODO: Add option to not compute summaries
  # TODO: Combine any gibbsm objects in ... with list
  # Interpret the which argument
  which <- match.arg(which)

  # Allow for either sequence of fits or list of fits, convert both to list
  fits <- if(missing(list)) {
    base::list(...)
  } else {
    list
  }

  # If user did not give names to the fits, add some default names
  default_names <- paste0("Fit ", seq_len(length(fits)))
  if(is.null(names(fits))) {
    names(fits) <- default_names
  }
  names(fits) <- ifelse(names(fits) == "", default_names, names(fits))

  if(!missing(summ)) { # Make sure summaries and fits are compatible
    if(!is(summ, "list")) {
      summ <- base::list(summ)
    }
    stopifnot(length(fits) == length(summ))
  } else { # Construct the summaries
    summ <- lapply(fits, function(f) summary(f))
  }

  # Below, we translate the coefficient parameter into a function "access" that allows us to access it in the fit object
  # In addition, is_interaction tells us whether or not the coefficient is an interaction coefficient (alpha/gamma) or not (beta)
  is_interaction <- TRUE
  if("gamma" == coefficient[1] & length(coefficient) == 1) { # The requested coefficient is gamma
    access <- function(obj) obj[["gamma"]]
  } else {
    is_alpha <- gsub("^alpha([0-9]*)", "\\1", coefficient)
    if("" == is_alpha[1] & length(is_alpha) == 1) { # it starts with alpha, but has no integer afterwards
      is_alpha <- 1 # by default, take first alpha
    }

    if(length(coefficient) == 1 & is_alpha[1] != coefficient[1]) { # it starts with alpha
      access <- function(obj) obj[["alpha"]][[as.numeric(is_alpha)]]
    } else {
      is_interaction <- FALSE
      is_beta <- gsub("^beta([0-9]*)", "\\1", coefficient)
      if("" == is_beta[1] & length(is_beta) == 1) { # it starts with beta, but has no integer afterwards
        # In this case, we want to get all environmental covariates involved, and use all
        list_covariates <- lapply(fits, function(fit) {
          colnames(as.matrix(fit$coefficients[["beta"]]))
        })
        is_beta <- as.character(unique(Reduce(c, list_covariates)))
      } else if(!all(is.na(suppressWarnings(as.numeric(is_beta))))) { # is_beta is an integer, telling us which beta to access
        is_beta <- as.numeric(is_beta)
      }

      # If it starts with beta, then this works. If not, assume that coefficient is an env covariate and access in this way too.
      access <- function(obj) {
        beta <- as.matrix(obj[["beta"]])
        if(is.character(is_beta)) {
          i <- intersect(is_beta, colnames(beta))
          z <- as.matrix(beta[, i])
          colnames(z) <- i
          z
        } else {
          beta
        }
      }
    }
  }


  # Extract list of alpha coefficients, each one corresponding to one of the fits
  estimates <- lapply(fits, function(f) {
    tryCatch(access(f$coefficients), error = function(err) NA)
  })

  # The fits could perhaps relate to fits which do not involve the same species. Construct a set of all species that appear in some of the fits
  species <- lapply(estimates, function(e) rownames(e))
  species <- unique(Reduce(c, species))

  # Species involving the focal ones
  species_subset <- if(!missing(involving)) {
    intersect(species, involving)
  } else {
    species
  }

  if(is_interaction) { # Create dataframe of two columns with all possible pairs of species
    df <- as.data.frame(expand.grid(from = species, to = species_subset, stringsAsFactors = FALSE))
    df <- df[!duplicated(t(apply(df, 1, sort))), ]

    if(which == "within") { # in this case, remove columns that are the same, so shows only between interactions
      df <- df[df$to == df$from, ]
    } else if(which == "between") {
      df <- df[df$to != df$from, ]
    }
  } else {
    df <- Reduce(rbind, lapply(is_beta, function(i) data.frame(from = species_subset, to = i)))
  }

  # At this point, df is a single dataframe containing all relevant pairs of species.
  # Since our estimates are a list of fits, convert df to a list of dataframe, each one corresponding to one of the fits.
  dfs <- lapply(seq_len(length(fits)), function(k) {
    d <- df
    d$lo <- sapply(seq_len(nrow(d)), function(i) { # Get the lower-endpoint of the CIs, apply function to sequence of numbers 1 to nrow.df
      tryCatch(access(summ[[k]]$lo)[d$from[i], d$to[i]], error = function(err) NA) # retrieve numbers from summ and apply to the rows
    })

    d$hi <- sapply(seq_len(nrow(d)), function(i) { # Get the upper-endpoint of the CIs
      tryCatch(access(summ[[k]]$hi)[d$from[i], d$to[i]], error = function(err) NA)
    })

    if(only_statistically_significant) { # in this case, remove non-stat significant
      d <- d[d$lo > 0 | d$hi < 0, ] # If non-statistically significant, remove row,. i.e. if low is > 0 or if high < 0, the row is kept
    }

    d$E <- sapply(seq_len(nrow(d)), function(i) { #these and below parts specify the values of point and bars
      tryCatch(estimates[[k]][d$from[i], d$to[i]], error = function(err) NA)
    })

    d$lo_numerical <- sapply(seq_len(nrow(d)), function(i) { # Lower-endpoint of the numerical CI (relating to the numerical error due to dummy points)
      tryCatch(access(summ[[k]]$lo_numerical)[d$from[i], d$to[i]], error = function(err) NA)
    })

    d$hi_numerical <- sapply(seq_len(nrow(d)), function(i) { # Upper-endpoint of the numerical CI (relating to the numerical error due to dummy points)
      tryCatch(access(summ[[k]]$hi_numerical)[d$from[i], d$to[i]], error = function(err) NA)
    })

    if(!is.null(full_names)) { # specifying names
      full_names <- as.list(full_names)

      d$from <- sapply(as.character(d$from), function(sp) if(length(full_names[[sp]]) == 0) sp else full_names[[sp]])
      if(is_interaction) {
        d$to <- sapply(as.character(d$to), function(sp) if(length(full_names[[sp]]) == 0) sp else full_names[[sp]])
      }
    }

    # Set name of the fit
    d$Fit <- names(fits)[k]

    d
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
  if(is_interaction) {
    levels(x) <- ifelse(df$from < df$to, paste0(df$from, enc2utf8(" \u2194 "), df$to), paste0(df$to, enc2utf8(" \u2194 "), df$from)) # Order lhs and rhs of interaction alphabetically
  } else { # In this, we are plotting env. covariates
    # If there is only one of them to plot, legend should just be the species
    # If not, we will be plotting multiple env. cov x species coefficients.
    if(length(is_beta) > 1) {
      levels(x) <- df$from
    } else {
      levels(x) <- df$from
    }
  }

  df$x <- factor(x, levels = levels(x)[order(df$average_estimates)])
  do_split <- FALSE
  if(!is_interaction) { # In this case, alphabetical order seems to be reasonable
    if(length(is_beta) > 1) {
      do_split <- TRUE
      df$split <- df$to
      df$x <- factor(x, levels = sort(levels(x), decreasing = TRUE))
    }
  }

  # Order Fit legend according to fits object
  df$Fit <- factor(df$Fit, levels = names(fits))

  # Remove NAs
  df <- df[rowSums(is.na(df)) == 0, ]

  # Set colours
  if(missing(colours)) {
    colours <- rep(c("black", "#5BBCD6", "#F2AD00", "#00A08A", "#FF0000"), length.out = length(fits))
  } else {
    colours <- rep(colours, length.out = length(fits))
  }

  g <- ggplot(df, aes_string(y = "x", x = "E")) # Visualisation

  if(highlight_zero) { # If this option is chosen, draw a red line before doing anything else
    g <- g + geom_vline(xintercept = 0, color = "red", linewidth = 1.5) # Vertical line at 0
  }

  g <- g +
    geom_point(aes_string(colour = "Fit"), size = 3, position = position_dodge(width = 0.2)) + # Size of point estimates
    geom_errorbar(aes_string(colour = "Fit",
                      xmin = "lo_numerical",
                      xmax = "hi_numerical"), width = 0, linewidth = 2, position = position_dodge(width = 0.2)) +
    geom_errorbar(aes_string(colour = "Fit",
                      xmin = "lo",
                      xmax = "hi"), width = 0.2, linewidth = 0.5, position = position_dodge(width = 0.2)) +
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

  # If xmin and xmax are given use those instead of defaults
  if(!missing(xmin) && !missing(xmax)) {
    g <- g + xlim(xmin, xmax) # Fix the x-axis range
  }

  # If there is a single fit, do not show legend
  if(length(fits) == 1) {
    g <- g + guides(colour = "none")
  }
  if(do_split) {
    g <- g + facet_grid(split ~ .)
  }
  g
}
#' Plot a `ppjsdm::gibbsm` object.
#'
#' This calls `ppjsdm::box_plot`, see the documentation for that function for details.
#' @method plot gibbsm
#' @param x A fit object obtained by calling `ppjsdm::gibbsm`.
#' @param ... Forwarded to `ppjsdm::box_plot`.
#' @export
#' @md
plot.gibbsm <- function(x, ...) {
  box_plot(x, ...)
}
