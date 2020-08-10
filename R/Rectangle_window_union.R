#' Rectangle window union constructor
#'
#' @param x_ranges Range along the x-axis.
#' @param y_ranges Range along the y-axis.
#' @importFrom methods is
#' @export
Rectangle_window_union <- local({
  function(x_ranges = list(c(0, 1)), y_ranges = list(c(0, 1))) {
    # Copy constructor
    if(nargs() == 1 && !missing(x_ranges) && is(x_ranges, "Rectangle_window_union")) {
      x_ranges
    } else { # Other constructors
      is_range <- function(r) {
        is.vector(r) &
        is.numeric(r) &
        length(r) == 2 &
        r[2L] >= r[1L]
      }
      if(!is.list(x_ranges) | !all(sapply(x_ranges, is_range))) {
        stop("x_ranges should be a list of numeric vectors of length 2 representing intervals (x_min, x_max).")
      }

      if(!is.list(y_ranges) | !all(sapply(y_ranges, is_range))) {
        stop("y_ranges should be a list of numeric vectors of length 2 representing intervals (x_min, x_max).")
      }

      # TODO: Make sure lists are the same length
      # TODO: Make sure they don't overlap
      structure(list(x_ranges = x_ranges, y_ranges = y_ranges), class = c("Rectangle_window_union"))
    }
  }
})

#' @method x_range Rectangle_window_union
#' @export
x_range.Rectangle_window_union <- function(window) {
  c(min(sapply(window$x_ranges, function(r) r[1])), max(sapply(window$x_ranges, function(r) r[2])))
}

#' @method y_range Rectangle_window_union
#' @export
y_range.Rectangle_window_union <- function(window) {
  c(min(sapply(window$y_ranges, function(r) r[1])), max(sapply(window$y_ranges, function(r) r[2])))
}

#' Return the area of a union of rectangle windows.
#'
#' @param window The window.
#' @export
area.Rectangle_window_union <- function(window) {
  sum(sapply(seq_len(length(window$x_ranges)), function(i) {
    x <- window$x_ranges[[i]]
    y <- window$y_ranges[[i]]
    (x[2] - x[1]) * (y[2] - y[1])
  }))
}

#' Return the volume of a rectangle window.
#'
#' @param window The window.
#' @export
volume.Rectangle_window_union <- function(window) {
  area.Rectangle_window_union(window)
}