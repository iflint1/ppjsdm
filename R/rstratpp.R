#' Sample a stratified binomial point processes
#'
#' @param window Simulation window. Default is a Rectangle window [0, 1]^2.
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
  if(is.list(delta_x)) {
    delta_x <- unlist(delta_x)
  }
  if(is.list(delta_y)) {
    delta_y <- unlist(delta_y)
  }

  rstratpp_cpp(window,
               delta_x,
               delta_y,
               nsim,
               types,
               drop,
               mark_range)
}
