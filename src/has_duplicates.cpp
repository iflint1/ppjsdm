#include <Rcpp.h>

#include "configuration/configuration_manipulation.hpp"
#include "configuration/configuration_wrapper.hpp"
#include "point/point_manipulation.hpp"

#include <algorithm> // std::sort, std::unique
#include <vector> // std::vector

//' Check if a configuration contains duplicates.
//'
//' @param configuration Configuration.
//' @export
//' @useDynLib ppjsdm, .registration = TRUE
//' @import Rcpp
// [[Rcpp::export]]
bool has_duplicates(Rcpp::List configuration) {
  const ppjsdm::Configuration_wrapper wrapped_configuration(configuration);
  const R_xlen_t number_points(ppjsdm::size(wrapped_configuration));
  std::vector<ppjsdm::Marked_point> points(number_points);
  for(R_xlen_t i(0); i < number_points; ++i) {
    points[i] = wrapped_configuration[i];
  }
  std::sort(points.begin(), points.end());
  const auto unique_end(std::unique(points.begin(), points.end()));

  return unique_end != points.end();
}
