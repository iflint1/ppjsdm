
#' Sample a binomial point processes
#'
#' @param window Observation window. Default is a rectangle window \eqn{[0, 1]^2}. Preferably a `ppjsdm` Window, such as `ppjsdm::Rectangle_window`, but also accepts `spatstat` `im` or `owin` objects.
#' @param n A vector representing the number of points of each types of the multipoint binomial point processes.
#' Default is a vector of same size as types, filled with ones.
#' @param nsim Number of samples to generate. Default is 1.
#' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
#' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration.
#' Default is TRUE.
#' @param mark_range Range of additional marks to give to the points.
#' @export
#' @md
rbinomialpp <- function(window = Rectangle_window(),
                 n = NULL,
                 nsim = 1,
                 types = NULL,
                 drop = TRUE,
                 mark_range = c(1.0, 1.0)) {
  n <- unlist(n)
  types <- unlist(types)
  number_types <- get_number_types_and_check_conformance(default_number_types = 1, types = types, n)$number_types
  n <- construct_if_missing(x = n, def = 1, nrows = number_types, matrix = FALSE)
  types <- make_types(types = types, size = number_types, might_contain_name = n)
  rbinomialpp_cpp(window = as.Window(window),
                  n = n,
                  nsim = nsim,
                  types = types,
                  drop = drop,
                  mark_range = mark_range)
}
