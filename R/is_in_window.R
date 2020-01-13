#' Check if point is in window
#'
#' @param x X-coordinates of points.
#' @param y Y-coordinates of points.
#' @param window Window.
#' @export
is_in_window <- function(x, y, window) UseMethod("is_in_window", window)

#' Check if points are in rectangle window
#'
#' @param x X-coordinates of points.
#' @param y Y-coordinates of points.
#' @param window Window
#' @export
is_in_window.Rectangle_window <- function(x, y, window) {
  stopifnot(length(x) == length(y))
  range <- cbind(x_range(window), y_range(window))
  points <- cbind(x, y)
  all(points >= range[1,] & points <= range[2,])
}
