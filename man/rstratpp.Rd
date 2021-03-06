% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstratpp.R
\name{rstratpp}
\alias{rstratpp}
\title{Sample a stratified binomial point processes}
\usage{
rstratpp(
  window = Rectangle_window(),
  delta_x = NULL,
  delta_y = NULL,
  nsim = 1,
  types = NULL,
  drop = TRUE,
  mark_range = c(1, 1)
)
}
\arguments{
\item{window}{Simulation window. Default is a Rectangle window [0, 1]^2.}

\item{nsim}{Number of samples to generate. Default is 1.}

\item{types}{Types of the points. Default is a vector (type1, type2, ...) of same size as n.}

\item{drop}{If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration.
Default is TRUE.}

\item{mark_range}{Range of additional marks to give to the points.}

\item{nx}{A vector representing the distance between tiles along the x axis, for each type.
Default is a vector of same size as types, filled with ones.}

\item{ny}{A vector representing the distance between tiles along the y axis, for each type.
Default is to set it to nx.}
}
\description{
Sample a stratified binomial point processes
}
