#' Rectangle window.
#'
#' @param x_range Range along the x-axis.
#' @param y_range Range along the y-axis.
#' @importFrom methods is
#' @export
#' @examples
#' # Different ways to construct a [0, 1] x [1, 3] rectangle window
#'
#' window <- ppjsdm::Rectangle_window(x_range = c(0, 1), y_range = c(1, 3))
#' print(window)
#'
#' window <- ppjsdm::Rectangle_window(y_range = c(1, 3))
#' print(window)
#'
#' window <- ppjsdm::Rectangle_window(window)
#' print(window)
#'
#' window <- ppjsdm::Rectangle_window(y_range = c(3, 1))
#' print(window)
Rectangle_window <- local({
  function(x_range = c(0, 1), y_range = c(0, 1)) {
    # Copy constructor
    if(nargs() == 1 && !missing(x_range) && is(x_range, "Rectangle_window")) {
      x_range
    } else { # Other constructors
      if(!is.vector(x_range) || !is.numeric(x_range) || length(x_range) != 2) {
        stop("x_range should be a numeric vector of length 2 representing an interval (x_min, x_max).")
      }

      if(!is.vector(y_range) || !is.numeric(y_range) || length(y_range) != 2) {
        stop("y_range should be a numeric vector of length 2 representing an interval (y_min, y_max).")
      }

      structure(list(x_range = c(min(x_range), max(x_range)),
                     y_range = c(min(y_range), max(y_range))), class = c("Rectangle_window"))
    }
  }
})

#' @method print Rectangle_window
#' @export
print.Rectangle_window <- function(x, ...) {
  str <- paste0("A rectangle window [", x$x_range[1], ", ", x$x_range[2], "] x [", x$y_range[1], ", ", x$y_range[2], "].")
  cat(str)
}

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
#' @param w Window.
#' @importFrom spatstat.geom area
#' @exportS3Method spatstat.geom::area Rectangle_window
#' @examples
#' # Construct a window
#'
#' window <- ppjsdm::Rectangle_window(x_range = c(0, 1), y_range = c(1, 3))
#'
#' # Compute area of a rectangle window
#'
#' print(area(window))
area.Rectangle_window <- function(w) {
  x <- w$x_range
  y <- w$y_range
  (x[2] - x[1]) * (y[2] - y[1])
}

#' Return the volume of a rectangle window.
#'
#' @param x Window.
#' @importFrom spatstat.geom volume
#' @exportS3Method spatstat.geom::volume Rectangle_window
#' @examples
#' # Construct a window
#'
#' window <- ppjsdm::Rectangle_window(x_range = c(0, 1), y_range = c(1, 3))
#'
#' # Compute volume of a rectangle window
#'
#' print(volume(window))
volume.Rectangle_window <- function(x) {
  area.Rectangle_window(x)
}

#' Convert a rectangle window to an owin from the `spatstat` package.
#'
#' @param W Window.
#' @param ... Unused.
#' @param fatal Unused.
#' @importFrom spatstat.geom as.owin owin
#' @exportS3Method spatstat.geom::as.owin Rectangle_window
#' @md
#' @examples
#' # Construct a window
#'
#' window <- ppjsdm::Rectangle_window(x_range = c(0, 1), y_range = c(1, 3))
#'
#' # Convert it to a spatstat object
#'
#' print(as.owin(window))
as.owin.Rectangle_window <- function(W, ..., fatal = TRUE) {
  owin(W$x_range, W$y_range)
}

