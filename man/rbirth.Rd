% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbirth.R
\name{rbirth}
\alias{rbirth}
\title{Sample a birth-death point process which has as its stationary distribution a multivariate saturated pairwise interaction Gibbs point process.
This function has been tested well but its theoretical foundation is still work in progress. The parameters are also temporary and might have different names in later versions.}
\usage{
rbirth(
  window,
  nsim = 1,
  types,
  horizon = 1,
  drop = TRUE,
  mark_range = c(1, 1),
  starting_configuration = NULL,
  nquad = 1000,
  nthreads = 1,
  birth_alpha,
  birth_gamma,
  birth_beta0,
  birth_covariates,
  birth_beta,
  birth_short_range,
  birth_medium_range,
  birth_long_range,
  birth_saturation,
  birth_model,
  birth_medium_range_model,
  death_alpha,
  death_gamma,
  death_beta0,
  death_covariates,
  death_beta,
  death_short_range,
  death_medium_range,
  death_long_range,
  death_saturation,
  death_model,
  death_medium_range_model
)
}
\arguments{
\item{window}{Observation window. Preferably a \code{ppjsdm} Window, such as \code{ppjsdm::Rectangle_window}, but also accepts \code{spatstat} \code{im} or \code{owin} objects.}

\item{nsim}{Number of samples to generate. Default is 1.}

\item{types}{Types of the points. Default is a vector (type1, type2, ...) of same size as n.}

\item{horizon}{Temporal horizon up to which we simulate}

\item{drop}{If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration. Default is TRUE.}

\item{mark_range}{Range of additional marks to give to the points.}

\item{starting_configuration}{Optional configuration to start with when using the birth-death process.}

\item{nquad}{Number of quadrature points used to approximate the integral of the birth rate.}

\item{nthreads}{Number of threads to use. Default is 1.}

\item{birth_alpha}{Repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.}

\item{birth_gamma}{Medium range repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.}

\item{birth_beta0}{A vector representing the log_intensities of the point processes.
Default is a vector of same size as types, filled with zeros}

\item{birth_covariates}{Covariates, with an empty list as a default.}

\item{birth_beta}{Fitted coefficients related to covariates. Default is square matrix of zeroes of the same
number of rows/columns as the covariates.}

\item{birth_short_range}{Symmetric matrix of short range interaction radii. Filled with 0.1 by default.}

\item{birth_medium_range}{Symmetric matrix of medium range interaction radii. Filled with 0 by default.}

\item{birth_long_range}{Symmetric matrix of long range interaction radii. Filled with 0 by default.}

\item{birth_saturation}{Saturation parameter of the point process. Default is 2.}

\item{birth_model}{String representing the model to use You can check the currently authorised models with a call to \code{show_short_range_models()}.}

\item{birth_medium_range_model}{String representing the model to use You can check the currently authorised models with a call to \code{show_medium_range_models()}.}

\item{death_alpha}{Repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.}

\item{death_gamma}{Medium range repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.}

\item{death_beta0}{A vector representing the log_intensities of the point processes.
Default is a vector of same size as types, filled with zeros}

\item{death_covariates}{Covariates, with an empty list as a default.}

\item{death_beta}{Fitted coefficients related to covariates. Default is square matrix of zeroes of the same
number of rows/columns as the covariates.}

\item{death_short_range}{Symmetric matrix of short range interaction radii. Filled with 0.1 by default.}

\item{death_medium_range}{Symmetric matrix of medium range interaction radii. Filled with 0 by default.}

\item{death_long_range}{Symmetric matrix of long range interaction radii. Filled with 0 by default.}

\item{death_saturation}{Saturation parameter of the point process. Default is 2.}

\item{death_model}{String representing the model to use You can check the currently authorised models with a call to \code{show_short_range_models()}.}

\item{death_medium_range_model}{String representing the model to use You can check the currently authorised models with a call to \code{show_medium_range_models()}.}
}
\description{
Sample a birth-death point process which has as its stationary distribution a multivariate saturated pairwise interaction Gibbs point process.
This function has been tested well but its theoretical foundation is still work in progress. The parameters are also temporary and might have different names in later versions.
}
