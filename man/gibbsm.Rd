% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gibbsm.R
\name{gibbsm}
\alias{gibbsm}
\title{Fit a multivariate Gibbs model to a dataset.}
\usage{
gibbsm(
  configuration_list,
  window,
  covariates,
  model,
  medium_range_model,
  short_range,
  medium_range,
  long_range,
  saturation,
  dummy = NULL,
  debug = FALSE,
  fitting_package = c("glm", "glmnet"),
  which_lambda = "AIC",
  max_dummy = 0,
  min_dummy = 0,
  dummy_factor = 0,
  nthreads = 1,
  use_regularization = FALSE,
  return_raw_fitting_data = FALSE,
  refit_glmnet = 0,
  dummy_distribution = "stratified",
  ...
)
}
\arguments{
\item{configuration_list}{A single configuration or a list of configurations assumed to be drawn from the multivariate Gibbs.}

\item{window}{Observation window. Preferably a \code{ppjsdm} Window, such as \code{ppjsdm::Rectangle_window}, but also accepts \code{spatstat} \code{im} or \code{owin} objects.}

\item{covariates}{Environmental covariates driving the intensity.}

\item{model}{String to represent the model we're calibrating. You can check the currently authorised models with a call to show_short_range_models().}

\item{medium_range_model}{String to represent the model we're calibrating. You can check the currently authorised models with a call to show_medium_range_models().}

\item{short_range}{Short range interaction radius. Filled with 0.1 by default.}

\item{medium_range}{Medium range interaction radius. Filled with 0 by default.}

\item{long_range}{Long range interaction radius. Filled with 0 by default.}

\item{saturation}{Saturation parameter of the point process.}

\item{dummy}{(Optional) dummy point configuration used in the fitting procedure.}

\item{debug}{Print debugging information?}

\item{fitting_package}{Which fitting package to use to solve the logistic regression?}

\item{which_lambda}{Which lambda to choose in the regularised fit? Choices are "AIC" (lambda that minimises AIC), "BIC" (lambda that minimises BIC) or "smallest" (smallest lambda).}

\item{max_dummy}{Maximum number of dummy points for each type. By default, follows the recommendation of Baddeley et al.}

\item{min_dummy}{Minimum number of dummy points for each type. By default, follows the recommendation of Baddeley et al.}

\item{dummy_factor}{How many more dummy points there are, compared to the points in the configuration. By default, follows the recommendation of Baddeley et al.}

\item{nthreads}{Maximum number of threads for parallel computing.}

\item{use_regularization}{Use the fitting package without regularization.}

\item{return_raw_fitting_data}{Return the raw fitting data, before calling the GLM fitting package. Mostly used for debugging purposes on large datasets.}

\item{refit_glmnet}{How many times more lambdas to use on a re-run of glmnet (currently only used with \code{fitting_package = glmnet}). For example, \code{refit_glmnet = 0.5} re-runs the fitting procedure with a 50\% longer lambda sequence.}

\item{dummy_distribution}{How should the dummy distribution be drawn? Choices are "binomial" or "stratified".}

\item{...}{Forwarded to the fitting package.}
}
\description{
Fit a multivariate Gibbs model to a dataset.
}
