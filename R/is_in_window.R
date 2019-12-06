# TODO: The more standard way is probably to have 2 columns and N rows, synch this with sample_from_window

#' Check if point is in window
#'
#' @param points Points.
#' @param window Window.
#' @export
setGeneric("is_in_window", function(points, window) {
  standardGeneric("is_in_window")
})

#' Check if points are in rectangle window
#'
#' @param points Points.
#' @param window Window
#' @export
setMethod("is_in_window", signature(window = "Rectangle_window"), function(points, window) {
  if(NROW(points) != 2) {
    stop("The points should have 2 rows and N columns, where N is the number of points and the rows are their coordinates.")
  }
  range <- cbind(x_range(window), y_range(window))
  all(points >= range[1,] & points <= range[2,])
})
