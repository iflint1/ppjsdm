% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chord_diagram.R
\name{chord_diagram_plot}
\alias{chord_diagram_plot}
\title{Plot a chord diagram of the coefficients of a fit object.}
\usage{
chord_diagram_plot(
  fit,
  coefficient = "alpha1",
  summ,
  only_statistically_significant = FALSE,
  full_names = NULL,
  compute_confidence_intervals = TRUE,
  classes = NULL,
  involving = NULL,
  how = c("only", "one"),
  include_self = TRUE,
  show_legend = TRUE,
  cex = 1,
  legend_cex = cex,
  ninteractions = 100,
  track_height = 0.05,
  show_bottom_left = TRUE,
  outward_facing_names = FALSE,
  show_grid_ticks = TRUE,
  sort_interactions = TRUE,
  only_between = FALSE,
  circle_margin = 0.1,
  big_gap = 5,
  repulsion_attraction_colours = c("blue", "red"),
  classes_colours
)
}
\arguments{
\item{fit}{A fit object obtained by a call to \code{ppjsdm::gibbsm}.}

\item{coefficient}{A string representing the coefficient to plot.
Choice of \code{alpha1}, \code{alpha2}, ..., to show one of the short-range interaction coefficients;
\code{alpha} to show all of the short-range potentials on the same chord diagram;
\code{gamma} to show the medium-range interaction coefficient.}

\item{summ}{Optional summary corresponding to the fit; if not provided it is obtained by calling \code{ppjsdm::summary.gibbsm}.}

\item{only_statistically_significant}{Only show statistically significant coefficients?}

\item{full_names}{Optional list of full names of types, if for example abbreviations were used when running the fit.}

\item{compute_confidence_intervals}{Compute the confidence intervals (which is slower) to highlight statistically significant interactions?}

\item{classes}{If this parameter is supplied, then colours are used to distinguish classes.
Should be a named vector/list, with the names corresponding to types, and the value equal to the class name.}

\item{involving}{Optional vector/list of types. Only coefficients involving these types will be plotted.}

\item{how}{If the \code{involving} argument is supplied, should it involve \emph{only} those types, or at least \emph{one} of those types (which is relevant if inter-type interactions are involved).}

\item{include_self}{Include self-interactions?}

\item{show_legend}{Show legend(s)?}

\item{cex}{Cex.}

\item{legend_cex}{Cex of the legend.}

\item{ninteractions}{Maximum number of interactions to include.}

\item{track_height}{Proportion of the chord diagram that the outer rim should occupy.}

\item{show_bottom_left}{Show the legend in the bottom-left? Only used if \code{classes} parameter is provided.}

\item{outward_facing_names}{Should the names of the types be outward facing or along the circle?}

\item{show_grid_ticks}{Show grid ticks on each of the sectors?}

\item{sort_interactions}{Should the interactions originating from a given type be sorted?}

\item{only_between}{Should we only plot interactions between different \emph{classes}? Only useful if classes were specified.}

\item{circle_margin}{Margin on the four sides of the circle.}

\item{big_gap}{Gap size between classes.}

\item{repulsion_attraction_colours}{Colours to represent repulsion and attraction.}

\item{classes_colours}{Colours to represent the different classes.}
}
\description{
The coefficients are shown in a chord diagram.
Highlighted opaque interactions are statistically significant.
The types are ordered according to their mean interaction, or alternatively sorted into classes if the \code{classes} argument is set.
}
\examples{
set.seed(1)

# Draw a configuration
configuration <- ppjsdm::rppp(lambda = c(A = 100, B = 100, C = 100, D = 100))

# Fit the data
fit <- ppjsdm::gibbsm(configuration)

# Chord diagram plot
ppjsdm::chord_diagram_plot(fit)
}
