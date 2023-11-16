#' Sample a stratified binomial point processes
#'
#' @param window Observation window. Default is a Rectangle window \eqn{[0, 1]^2}. Preferably a `ppjsdm` Window, such as `ppjsdm::Rectangle_window`, but also accepts `spatstat` `im` or `owin` objects.
#' @param delta_x A vector representing the distance between tiles along the x axis, for each type.
#' Default is a vector of same size as types, filled with ones.
#' @param delta_y A vector representing the distance between tiles along the y axis, for each type.
#' Default is to set it to delta_x
#' @param nsim Number of samples to generate. Default is 1.
#' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
#' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration.
#' Default is TRUE.
#' @param mark_range Range of additional marks to give to the points.
#' @export
rstratpp <- function(window = Rectangle_window(),
                     delta_x = NULL,
                     delta_y = NULL,
                     nsim = 1,
                     types = NULL,
                     drop = TRUE,
                     mark_range = c(1.0, 1.0)) {
  delta_x <- unlist(delta_x)
  delta_y <- unlist(delta_y)
  types <- unlist(types)
  number_types <- get_number_types_and_check_conformance(default_number_types = 1, types = types, delta_x, delta_y)$number_types
  delta_x <- construct_if_missing(x = delta_x, def = 1, nrows = number_types, matrix = FALSE)
  if(is.null(delta_y)) {
    delta_y <- delta_x
  }
  types <- make_types2(types = types, size = number_types, delta_x, delta_y)

  rstratpp_cpp(window = as.Window(window),
               delta_x = unlist(delta_x),
               delta_y = unlist(delta_y),
               nsim = nsim,
               types = types,
               drop = drop,
               mark_range = mark_range)
}
