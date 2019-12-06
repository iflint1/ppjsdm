#' A disk window
#'
#' @param centre Centre of the disk.
#' @param radius Radius of the disk.
#' @export
setClass("Disk_window", slots = list(centre = "numeric", radius = "numeric"))

#' Disk window constructor
#'
#' @param centre Centre of the disk.
#' @param radius Radius of the disk.
#' @importFrom methods is new
#' @export
Disk_window <- local({

  function(centre = c(0, 0), radius = 1) {
    # Copy constructor
    if(nargs() == 1 && !missing(centre) && is(centre, "Disk_window")) {
      centre
    } else { # Other constructors
      if(!is.vector(centre) || length(centre) != 2) {
        stop("centre should be a vector of length 2 representing the coordinates of the centre.")
      }
      if(!is.vector(radius) || length(radius) != 1 || radius < 0) {
        stop("radius should be a non-negative vector of length 1 representing the radius.")
      }

      new("Disk_window", centre = centre, radius = radius)
    }
  }

})

#' Access centre of a disk window
#'
#' @param window The window
#' @export
centre <- function(window) {
  window@centre
}

#' Access radius of a disk window
#'
#' @param window The window
#' @export
radius <- function(window) {
  window@radius
}
