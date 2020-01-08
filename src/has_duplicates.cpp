#include <Rcpp.h>

#include "configuration_utilities.h"

#include <algorithm> // std::sort, std::unique
#include <tuple> // std::tuple
#include <vector> // std::vector

using Marked_point = std::tuple<double, double, int>;

//' Check if a configuration contains duplicates.
//'
//' @param configuration The configuration.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
bool has_duplicates(Rcpp::List configuration) {
  const Configuration_wrapper wrapped_configuration(configuration);
  const R_xlen_t number_points(wrapped_configuration.get_number_points());
  std::vector<Marked_point> points(number_points);
  for(R_xlen_t i(0); i < number_points; ++i) {
    points[i] = Marked_point(wrapped_configuration.x()[i], wrapped_configuration.y()[i], wrapped_configuration.types()[i]);
  }
  std::sort(points.begin(), points.end());
  const auto unique_end{std::unique(points.begin(), points.end())};

  return unique_end != points.end();
}
