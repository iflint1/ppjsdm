#include <Rcpp.h>

#include "utility/configuration_manipulation.h"
#include "utility/configuration_wrappers.h"

#include <algorithm> // std::sort, std::unique
#include <tuple> // std::tuple
#include <vector> // std::vector

namespace ppjsdm {
namespace detail {

using Marked_point = std::tuple<double, double, int>;

} // namespace detail
} // namespace ppjsdm

//' Check if a configuration contains duplicates.
//'
//' @param configuration Configuration.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
bool has_duplicates(Rcpp::List configuration) {
  using Marked_point = ppjsdm::detail::Marked_point;
  const ppjsdm::Configuration_wrapper wrapped_configuration(configuration);
  const R_xlen_t number_points(ppjsdm::size(wrapped_configuration));
  std::vector<Marked_point> points(number_points);
  for(R_xlen_t i(0); i < number_points; ++i) {
    points[i] = Marked_point(wrapped_configuration.x(i), wrapped_configuration.y(i), wrapped_configuration.types(i));
  }
  std::sort(points.begin(), points.end());
  const auto unique_end(std::unique(points.begin(), points.end()));

  return unique_end != points.end();
}
