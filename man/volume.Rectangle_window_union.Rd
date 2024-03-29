% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rectangle_window_union.R
\name{volume.Rectangle_window_union}
\alias{volume.Rectangle_window_union}
\title{Return the volume of a rectangle window.}
\usage{
\method{volume}{Rectangle_window_union}(x)
}
\arguments{
\item{x}{Window.}
}
\description{
Return the volume of a rectangle window.
}
\examples{
# Construct a window

window <- ppjsdm::Rectangle_window_union(x_ranges = list(c(0, 1), c(1, 2)), y_ranges = list(c(1, 3), c(2, 3)))

# Compute volume of a union of rectangle windows

print(volume(window))
}
