#' A rectangle window
#'
#' @param x_range Range along the x-axis.
#' @param y_range Range along the y-axis.
#' @export
setClass("Rectangle_window", slots = list(x_range = "numeric", y_range = "numeric"))

#' Rectangle window constructor
#'
#' @param x_range Range along the x-axis.
#' @param y_range Range along the y-axis.
#' @importFrom methods is new
#' @export
Rectangle_window <- local({

  function(x_range = c(0, 1), y_range = c(0, 1)) {
    # Copy constructor
    if(nargs() == 1 && !missing(x_range) && is(x_range, "Rectangle_window")) {
      x_range
    } else { # Other constructors
      if(!is.vector(x_range) || length(x_range) != 2 || x_range[2L] < x_range[1L]) {
        stop("x_range should be a vector of length 2 representing an interval (x_min, x_max).")
      }

      if(!is.vector(y_range) || length(y_range) != 2 || y_range[2L] < y_range[1L]) {
        stop("y_range should be a vector of length 2 representing an interval (y_min, y_max).")
      }

      new("Rectangle_window", x_range = x_range, y_range = y_range)
    }
  }

})

#' Access range along the x-axis of a rectangle window
#'
#' @param window The window.
#' @export
x_range <- function(window) {
  window@x_range
}

#' Access range along the y-axis of a rectangle window
#'
#' @param window The window.
#' @export
y_range <- function(window) {
  window@y_range
}

#' Return the volume of the window.
#'
#' @param window The window.
#' @export
setMethod("window_volume", signature(window = "Rectangle_window"), function(window) {
  x <- window@x_range
  y <- window@y_range
  (x[2] - x[1]) * (y[2] - y[1])
})

#' Convert a rectangle window to an owin from the SpatStat package.
#'
#' @param W Window.
#' @param ... Unused.
#' @param fatal Unused.
#' @importFrom spatstat owin
#' @export
as.owin.Rectangle_window <- function(W, ..., fatal = TRUE) {
  owin(W@x_range, W@y_range)
}

