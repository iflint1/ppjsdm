
#' Sample a Poisson point processes
#'
#' @param window Observation window. Default is a Rectangle window \eqn{[0, 1]^2}. Preferably a `ppjsdm` Window, such as `ppjsdm::Rectangle_window`, but also accepts `spatstat` `im` or `owin` objects.
#' @param lambda A vector representing the intensities of the multipoint Poisson point processes.
#' Default is a vector of same size as types, filled with 100s.
#' @param nsim Number of samples to generate. Default is 1.
#' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
#' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration.
#' Default is TRUE.
#' @param mark_range Range of additional marks to give to the points.
#' @export
#' @md
rppp <- function(window = Rectangle_window(),
                 lambda = NULL,
                 nsim = 1,
                 types = NULL,
                 drop = TRUE,
                 mark_range = c(1.0, 1.0)) {
  lambda <- unlist(lambda)
  types <- unlist(types)
  number_types <- get_number_types_and_check_conformance(default_number_types = 1, lambda, types)
  lambda <- construct_if_missing(x = lambda, def = 100, nrows = number_types, matrix = FALSE)
  types <- make_types(types = types, size = number_types, might_contain_name = lambda)

  rppp_cpp(window = as.Window(window),
           lambda = lambda,
           nsim = nsim,
           types = types,
           drop = drop,
           mark_range = mark_range)
}
