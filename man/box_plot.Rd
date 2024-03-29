% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/box_plot.R
\name{box_plot}
\alias{box_plot}
\title{Plot the coefficients of fit objects.}
\usage{
box_plot(
  ...,
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
  xmax
)
}
\arguments{
\item{...}{Any number of fit objects obtained by a call to \code{ppjsdm::gibbsm}.}

\item{list}{Some more fits provided as a list.}

\item{coefficient}{A string representing the coefficient to plot.
Choice of \code{alpha1}, \code{alpha2}, ..., to show one of the short-range interaction coefficients;
\code{gamma} to show the medium-range interaction coefficient;
\code{beta1}, \code{beta2}, ... or actual name of the covariate for one of the regression coefficients.
\code{alpha} shows a facet plot of all short-range interaction coefficients, while \code{beta} does the same for the regression coefficients.}

\item{summ}{Optional list of summaries corresponding to the fits; if not provided they are obtained by calling \code{ppjsdm::summary.gibbsm}}

\item{only_statistically_significant}{Only show statistically significant coefficients?}

\item{which}{If plotting interaction coefficients, which ones do we want to plot?}

\item{full_names}{Optional list of full names of types, if for example abbreviations were used when running the fit.}

\item{compute_confidence_intervals}{Compute the confidence intervals (which is slower) or just show the point estimates?}

\item{classes}{If this parameter is supplied, then colours are used to distinguish classes instead of fits.
Only works when given a single fit.
Should be a named vector/list, with the names corresponding to types, and the value equal to the class name.}

\item{involving}{Optional vector/list of types. Only coefficients involving these types will be plotted.}

\item{how}{If the \code{involving} argument is supplied, should it involve \emph{only} those types, or at least \emph{one} of those types (which is relevant if inter-type interactions are involved).}

\item{title}{Plot title.}

\item{colours}{Optional vector of colours to represent the different fits/classes.}

\item{highlight_zero}{Highlight the zero value with a red line?}

\item{text_size}{Text size.}

\item{base_size}{Base size.}

\item{xmin}{Optional plot minimum x.}

\item{xmax}{Optional plot maximum x.}
}
\description{
The coefficients and corresponding confidence intervals are shown in a box-like plot.
The inner thick parts of the error bars represent numerical uncertainty, due to the number and distribution of dummy points.
The outer thin part of the error bars represent the total statistical + numerical uncertainty.
Formally, these are theoretical asymptotic 95\% confidence intervals.
As a convenience, this function can alternatively be applied by calling \code{plot} on a fit object.
}
\examples{
set.seed(1)

# Construct a configuration

configuration <- ppjsdm::rppp(lambda = c(A = 500, B = 500))

# Fit the model

fit <- ppjsdm::gibbsm(configuration, covariates = list(x = function(x, y) x))

# Plot the coefficients

ppjsdm::box_plot(fit)

ppjsdm::box_plot(fit, coefficient = "x")

}
