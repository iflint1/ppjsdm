#' Im window.
#'
#' @param im `spatstat` `im` object.
#' @param window Optional `spatstat` window object.
#' @importFrom methods is
#' @importFrom spatstat.geom as.owin is.im
#' @export
#' @md
#' @examples
#' if(require(spatstat.geom)) { # If no spatstat, im windows cannot be constructed
#' # Construct an im window
#'
#' window <- ppjsdm::Im_window(as.im.function(function(x, y) ifelse(x^2 + y^2 <= 1, 1, NA), W = owin()))
#'
#' print(window)
#'
#' plot(window)
#' }
Im_window <- local({
  function(im, window) {
    # Copy constructor
    if(nargs() == 1 && !missing(im) && is(im, "Im_window")) {
      im
    } else { # Other constructors
      if(!is.im(im)) {
        im <- if(!missing(window)) {
          as.im(im, W = as.owin(window))
        } else {
          as.im(im)
        }
      }
      structure(list(im = im), class = c("Im_window", "Window"))
    }
  }
})

#' @importFrom spatstat.geom print.im
#' @method print Im_window
#' @export
print.Im_window <- function(x, ...) {
  cat("An im window: \n")
  print.im(x$im)
}

#' @method x_range Im_window
#' @export
x_range.Im_window <- function(window) {
  window$im$xrange
}

#' @method y_range Im_window
#' @export
y_range.Im_window <- function(window) {
  window$im$yrange
}

#' Return the area of an im window.
#'
#' @param w Window.
#' @importFrom spatstat.geom area area.owin
#' @exportS3Method spatstat.geom::area Im_window
#' @examples
#' if(require(spatstat.geom)) { # If no spatstat, im windows cannot be constructed
#' # Construct a window
#'
#' window <- ppjsdm::Im_window(as.im(function(x, y) ifelse(x^2 + y^2 <= 1, 1, NA), W = owin()))
#'
#' # Compute area of an im window
#'
#' print(area(window))
#' }
area.Im_window <- function(w) {
  area.owin(as.owin.Im_window(w))
}

#' Return the volume of an im window.
#'
#' @param x Window.
#' @importFrom spatstat.geom volume
#' @exportS3Method spatstat.geom::volume Im_window
#' @examples
#' if(require(spatstat.geom)) { # If no spatstat, im windows cannot be constructed
#' # Construct a window
#'
#' window <- ppjsdm::Im_window(as.im(function(x, y) ifelse(x^2 + y^2 <= 1, 1, NA), W = owin()))
#'
#' # Compute volume of an im window
#'
#' print(volume(window))
#' }
volume.Im_window <- function(x) {
  area.Im_window(x)
}

#' Convert an im window to an owin from the `spatstat` package.
#'
#' @param W Window.
#' @param ... Unused.
#' @param fatal Unused.
#' @importFrom spatstat.geom as.owin as.owin.im
#' @exportS3Method spatstat.geom::as.owin Im_window
#' @md
#' @examples
#' if(require(spatstat.geom)) { # If no spatstat, im windows cannot be constructed
#' # Construct a window
#'
#' window <- ppjsdm::Im_window(as.im(function(x, y) ifelse(x^2 + y^2 <= 1, 1, NA), W = owin()))
#'
#' # Convert it to a spatstat object
#'
#' print(as.owin(window))
#' }
as.owin.Im_window <- function(W, ..., fatal = TRUE) {
  as.owin.im(W$im)
}
