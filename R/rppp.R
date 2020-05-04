
#' Sample a Poisson point processes
#'
#' @param window Simulation window. Default is a Rectangle window [0, 1]^2.
#' @param lambda A vector representing the intensities of the multipoint Poisson point processes.
#' Default is a vector of same size as types, filled with ones.
#' @param nsim Number of samples to generate. Default is 1.
#' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
#' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration.
#' Default is TRUE.
#' @param mark_range Range of additional marks to give to the points.
#' @export
rppp <- function(window = Rectangle_window(),
                 lambda = NULL,
                 nsim = 1,
                 types = NULL,
                 drop = TRUE,
                 mark_range = c(1.0, 1.0)) {
  if(is.list(lambda)) {
    lambda <- unlist(lambda)
  }
  rppp_cpp(window,
           lambda,
           nsim,
           types,
           drop,
           mark_range)
}
