#' Rectangle window union.
#'
#' @param x_ranges List of rectangle ranges along the x-axis.
#' @param y_ranges List of rectangle ranges along the y-axis.
#' @importFrom methods is
#' @export
#' @examples
#' # Different ways to construct a union of a [0, 1] x [1, 3] rectangle window
#' # and a [1, 2] x [2, 3] rectangle window.
#'
#' window <- ppjsdm::Rectangle_window_union(x_ranges = list(c(0, 1), c(1, 2)), y_ranges = list(c(1, 3), c(2, 3)))
#' print(window)
#'
#' window <- ppjsdm::Rectangle_window_union(x_ranges = list(c(0, 1), c(1, 2)), y_ranges = list(c(3, 1), c(3, 2)))
#' print(window)
#'
#' window <- ppjsdm::Rectangle_window_union(window)
#' print(window)
Rectangle_window_union <- local({
  function(x_ranges = list(c(0, 1)), y_ranges = list(c(0, 1))) {
    # Copy constructor
    if(nargs() == 1 && !missing(x_ranges) && is(x_ranges, "Rectangle_window_union")) {
      x_ranges
    } else { # Other constructors
      is_range <- function(r) {
        is.vector(r) &
        is.numeric(r) &
        length(r) == 2
      }
      if(!is.list(x_ranges) | !all(sapply(x_ranges, is_range))) {
        stop("x_ranges should be a list of numeric vectors of length 2 representing intervals (x_min, x_max).")
      }

      if(!is.list(y_ranges) | !all(sapply(y_ranges, is_range))) {
        stop("y_ranges should be a list of numeric vectors of length 2 representing intervals (x_min, x_max).")
      }

      if(length(x_ranges) != length(y_ranges)) {
        stop("The list of x/y ranges should be of the same length.")
      }

      # make sure none of the ranges have inverted intervals
      x_ranges <- lapply(x_ranges, function(x) c(min(x), max(x)))
      y_ranges <- lapply(y_ranges, function(x) c(min(x), max(x)))

      # Make sure none of the rectangles overlap
      if(length(x_ranges) > 1) {
        for(i in 1:(length(x_ranges) - 1)) {
          for(j in (i + 1):length(x_ranges)) {
            if(x_ranges[[i]][1] < x_ranges[[j]][2] &&
               x_ranges[[i]][2] > x_ranges[[j]][1] &&
               y_ranges[[i]][2] > y_ranges[[j]][1] &&
               y_ranges[[i]][1] < y_ranges[[j]][2]) {
              stop(paste0("Two of the given rectangles overlap, first with xs (", paste0(x_ranges[[i]], collapse = ", "),
              ") and ys (", paste0(y_ranges[[i]], collapse = ", "), "), second with first with xs (",
              paste0(x_ranges[[j]], collapse = ", "), ") and ys (", paste0(y_ranges[[j]], collapse = ", "), ")."))
            }
          }
        }
      }

      structure(list(x_ranges = x_ranges, y_ranges = y_ranges), class = c("Rectangle_window_union"))
    }
  }
})

#' @method print Rectangle_window_union
#' @export
print.Rectangle_window_union <- function(x, ...) {
  str <- "A union of rectangle windows: \n"
  for(i in seq_len(length(x$x_ranges))) {
    str <- paste0(str, "(", i, ") [", x$x_ranges[[i]][1], ", ", x$x_ranges[[i]][2], "] x [",
                  x$y_ranges[[i]][1], ", ", x$y_ranges[[i]][2], "]")
    if(i < length(x$x_ranges)) {
      str <- paste0(str, ",\n")
    } else {
      str <- paste0(str, ".")
    }
  }
  cat(str)
}

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
#' @param w Window.
#' @importFrom spatstat.geom area
#' @export
#' @examples
#' # Construct a window
#'
#' window <- ppjsdm::Rectangle_window_union(x_ranges = list(c(0, 1), c(1, 2)), y_ranges = list(c(1, 3), c(2, 3)))
#'
#' # Compute area of a union of rectangle windows
#'
#' print(area(window))
area.Rectangle_window_union <- function(w) {
  sum(sapply(seq_len(length(w$x_ranges)), function(i) {
    x <- w$x_ranges[[i]]
    y <- w$y_ranges[[i]]
    (x[2] - x[1]) * (y[2] - y[1])
  }))
}

#' Return the volume of a rectangle window.
#'
#' @param x Window.
#' @importFrom spatstat.geom volume
#' @export
#' @examples
#' # Construct a window
#'
#' window <- ppjsdm::Rectangle_window_union(x_ranges = list(c(0, 1), c(1, 2)), y_ranges = list(c(1, 3), c(2, 3)))
#'
#' # Compute volume of a union of rectangle windows
#'
#' print(volume(window))
volume.Rectangle_window_union <- function(x) {
  area.Rectangle_window_union(x)
}

#' Convert a union of rectangle windows to an owin from the `spatstat` package.
#'
#' @param W Window.
#' @param ... Unused.
#' @param fatal Unused.
#' @importFrom spatstat.geom owin union.owin
#' @export
#' @md
#' @examples
#' # Construct a window
#'
#' window <- ppjsdm::Rectangle_window_union(x_ranges = list(c(0, 1), c(1, 2)), y_ranges = list(c(1, 3), c(2, 3)))
#'
#' # Convert it to a spatstat object
#'
#' print(as.owin(window))
#' plot(as.owin(window))
as.owin.Rectangle_window_union <- function(W, ..., fatal = TRUE) {
  owins <- lapply(seq_len(length(W$x_ranges)), function(i) owin(W$x_ranges[[i]], W$y_ranges[[i]]))
  do.call(union.owin, owins)
}
