#' Range along x-axis of enclosing rectangle of a window.
#'
#' @param window Window.
#' @export
#' @examples
#' # Construct a disk window
#'
#' window <- ppjsdm::Disk_window(centre = c(2, 1), radius = 2)
#'
#' # Compute the range along the x-axis
#'
#' print(ppjsdm::x_range(window))
#'
x_range <- function(window) {
  UseMethod("x_range", window)
}

#' Range along y-axis of enclosing rectangle of a window.
#'
#' @param window Window.
#' @export
#' @examples
#' # Construct a disk window
#'
#' window <- ppjsdm::Disk_window(centre = c(2, 1), radius = 2)
#'
#' # Compute the range along the y-axis
#'
#' print(ppjsdm::y_range(window))
#'
y_range <- function(window) {
  UseMethod("y_range", window)
}
