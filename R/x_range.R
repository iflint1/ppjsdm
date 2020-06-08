#' Range along x-axis of enclosing rectangle of window.
#'
#' @param window Window.
#' @export
x_range <- function(window) {
  UseMethod("x_range", window)
}

#' Range along y-axis of enclosing rectangle of window.
#'
#' @param window Window.
#' @export
y_range <- function(window) {
  UseMethod("y_range", window)
}
