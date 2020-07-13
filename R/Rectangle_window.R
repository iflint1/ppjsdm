#' Rectangle window constructor
#'
#' @param x_range Range along the x-axis.
#' @param y_range Range along the y-axis.
#' @importFrom methods is
#' @export
Rectangle_window <- local({
  function(x_range = c(0, 1), y_range = c(0, 1)) {
    # Copy constructor
    if(nargs() == 1 && !missing(x_range) && is(x_range, "Rectangle_window")) {
      x_range
    } else { # Other constructors
      if(!is.vector(x_range) || !is.numeric(x_range) || length(x_range) != 2 || x_range[2L] < x_range[1L]) {
        stop("x_range should be a numeric vector of length 2 representing an interval (x_min, x_max).")
      }

      if(!is.vector(y_range) || !is.numeric(y_range) || length(y_range) != 2 || y_range[2L] < y_range[1L]) {
        stop("y_range should be a numeric vector of length 2 representing an interval (y_min, y_max).")
      }

      structure(list(x_range = x_range, y_range = y_range), class = c("Rectangle_window"))
    }
  }
})

#' @method x_range Rectangle_window
#' @export
x_range.Rectangle_window <- function(window) {
  window$x_range
}

#' @method y_range Rectangle_window
#' @export
y_range.Rectangle_window <- function(window) {
  window$y_range
}

#' Return the area of a rectangle window.
#'
#' @param window The window.
#' @export
area.Rectangle_window <- function(window) {
  x <- window$x_range
  y <- window$y_range
  (x[2] - x[1]) * (y[2] - y[1])
}

#' Return the volume of a rectangle window.
#'
#' @param window The window.
#' @export
volume.Rectangle_window <- function(window) {
  area.Rectangle_window(window)
}

#' Convert a rectangle window to an owin from the SpatStat package.
#'
#' @param W Window.
#' @param ... Unused.
#' @param fatal Unused.
#' @importFrom spatstat owin
#' @export
as.owin.Rectangle_window <- function(W, ..., fatal = TRUE) {
  owin(W$x_range, W$y_range)
}

