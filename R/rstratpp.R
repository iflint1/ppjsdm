#' Sample a stratified binomial point processes
#'
#' @param window Simulation window. Default is a Rectangle window [0, 1]^2.
#' @param nx A vector representing the number of tiles in each column, for each type.
#' Default is a vector of same size as types, filled with ones.
#' @param ny A vector representing the number of tiles in each row, for each type.
#' Default is to set it to nx.
#' @param nsim Number of samples to generate. Default is 1.
#' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
#' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration.
#' Default is TRUE.
#' @param mark_range Range of additional marks to give to the points.
#' @export
rstratpp <- function(window = Rectangle_window(),
                     nx = NULL,
                     ny = NULL,
                     nsim = 1,
                     types = NULL,
                     drop = TRUE,
                     mark_range = c(1.0, 1.0)) {
  if(is.list(nx)) {
    nx <- unlist(nx)
  }
  if(is.list(ny)) {
    ny <- unlist(ny)
  }

  rstratpp_cpp(window,
                  nx,
                  ny,
                  nsim,
                  types,
                  drop,
                  mark_range)
}
