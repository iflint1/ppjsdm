% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.gibbsm.R
\name{summary.gibbsm}
\alias{summary.gibbsm}
\title{Summary of a Fitted Gibbs point process.
The function can be supplied multiple fits, in which case the aggregate estimator is computed.}
\usage{
\method{summary}{gibbsm}(
  object,
  ...,
  list,
  debug = FALSE,
  time_limit = Inf,
  unit = "hours",
  nthreads = NULL,
  npoints = 2000,
  multiple_windows = TRUE,
  assume_fitted_to_same_data = FALSE
)
}
\arguments{
\item{object}{A fitted model object.}

\item{...}{Other fitted model objects.}

\item{list}{List of fits to consider, in addition to object.}

\item{debug}{Display debug information?}

\item{time_limit}{Time limit measured in  `unit` that can be spent running this function. This is a rough attempted upper-bound on run
time that is sometimes exceeded if the first batch already takes longer than this limit.
Nb: time_limit / 4 is allocated for the computation of each of the 4 underlying matrices that constitute the variance-covariance matrix.
In mid-sized datasets, most of the time is spent computing A2 + A3 and (to a lesser extent) G2. For such data, A1 and S may be
computed almost instantly, and time_limit / 2 will remain to compute the remaining two matrices. Use debug = TRUE to
assess this in the target use-case.}

\item{unit}{Unit used to measure the time limit (hours, mins, secs, etc).}

\item{nthreads}{(optional) number of threads to use.}

\item{npoints}{Target number of points in the restricted window that the vcov matrix is computed on.
Computation is slower for larger values, but the vcov matrix is then better approximated. The larger this is, the better.}

\item{multiple_windows}{Compute A2 and A3 on a lot of small windows and which are then averaged out, or only on a single restricted window?
Setting this to FALSE uses the standard estimate, but can be excessively slow to compute on large datasets.
The default of TRUE attempts to compute the estimate on multiple windows containing around npoints data points in each.}

\item{assume_fitted_to_same_data}{Should the function assume that the model was fitted to the same data?
If so and if multiple fits are supplied, the function assumes that all configurations are identical.
If not, the function assumes that the different fits were obtained through different realisations of the same point process.}
}
\description{
Summary of a Fitted Gibbs point process.
The function can be supplied multiple fits, in which case the aggregate estimator is computed.
}
