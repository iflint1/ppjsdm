#' Disk window constructor
#'
#' @param centre Centre of the disk.
#' @param radius Radius of the disk.
#' @importFrom methods is
#' @export
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

      structure(list(centre = centre, radius = radius), class = c("Disk_window", "Window"))
    }
  }
})

#' Access centre of a disk window
#'
#' @param window The window
#' @export
centre <- function(window) {
  window$centre
}

#' Access radius of a disk window
#'
#' @param window The window
#' @export
radius <- function(window) {
  window$radius
}

#' Return the volume of the window.
#'
#' @param window The window.
#' @export
window_volume <- function(window) UseMethod("window_volume", window)

#' Return the volume of the window.
#'
#' @param window The window.
#' @export
window_volume.Disk_window <- function(window) {
  r <- window$radius
  pi * r * r
}
