#' Disk window.
#'
#' @param centre Centre of the disk.
#' @param radius Radius of the disk.
#' @importFrom methods is
#' @export
#' @examples
#' # Construct a disk window with centre (0, 2) and radius 2.
#'
#' window <- ppjsdm::Disk_window(centre = c(0, 2), radius = 2)
#' print(window)
#'
#' window <- ppjsdm::Disk_window(window)
#' print(window)
Disk_window <- local({
  function(centre = c(0, 0), radius = 1) {
    # Copy constructor
    if(nargs() == 1 && !missing(centre) && is(centre, "Disk_window")) {
      centre
    } else { # Other constructors
      if(!is.vector(centre) || !is.numeric(centre) || length(centre) != 2) {
        stop("centre should be a numeric vector of length 2 representing the coordinates of the centre.")
      }
      if(!is.vector(radius) || !is.numeric(radius) || length(radius) != 1 || radius < 0) {
        stop("radius should be a non-negative numeric vector of length 1 representing the radius.")
      }

      structure(list(centre = centre, radius = radius), class = c("Disk_window"))
    }
  }
})

#' @method print Disk_window
#' @export
print.Disk_window <- function(x, ...) {
  str <- paste0("A disk window with centre (", x$centre[1], ", ", x$centre[2], ") and radius ", x$radius, ".")
  cat(str)
}

#' @method x_range Disk_window
#' @export
x_range.Disk_window <- function(window) {
  r <- window$radius
  x <- window$centre[1]
  c(x - r, x + r)
}

#' @method y_range Disk_window
#' @export
y_range.Disk_window <- function(window) {
  r <- window$radius
  y <- window$centre[2]
  c(y - r, y + r)
}

#' Return the area of a disk window.
#'
#' @param w Window.
#' @importFrom spatstat.geom area
#' @export
#' @examples
#' # Construct a window
#'
#' window <- ppjsdm::Disk_window(centre = c(0, 2), radius = 2)
#'
#' # Compute area of a rectangle window
#'
#' print(area(window))
area.Disk_window <- function(w) {
  r <- w$radius
  pi * r * r
}

#' Return the volume of a disk window.
#'
#' @param x Window.
#' @importFrom spatstat.geom volume
#' @export
#' @examples
#' # Construct a window
#'
#' window <- ppjsdm::Disk_window(centre = c(0, 2), radius = 2)
#'
#' # Compute volume of a rectangle window
#'
#' print(volume(window))
volume.Disk_window <- function(x) {
  area.Disk_window(x)
}

#' Convert a disk window to an owin from the SpatStat package.
#'
#' @param W Window.
#' @param ... Unused.
#' @param fatal Unused.
#' @importFrom spatstat.geom disc
#' @export
#' @examples
#' # Construct a window
#'
#' window <- ppjsdm::Disk_window(centre = c(0, 2), radius = 2)
#'
#' # Convert it to a spatstat object
#'
#' print(as.owin(window))
as.owin.Disk_window <- function(W, ..., fatal = TRUE) {
  disc(radius = W$radius, centre = W$centre)
}
