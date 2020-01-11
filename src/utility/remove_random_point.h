#ifndef INCLUDE_PPJSDM_REMOVE_RANDOM_POINT
#define INCLUDE_PPJSDM_REMOVE_RANDOM_POINT

#include "configuration_manipulation.h"

#include <iterator> // std::iterator_traits

namespace ppjsdm {

// Assumes that the configuration is non-empty
template<typename Configuration>
inline auto remove_random_point(Configuration& configuration) {
  using difference_type = typename std::iterator_traits<decltype(configuration.begin())>::difference_type;
  const difference_type index(Rcpp::sample(size(configuration), 1, false, R_NilValue, false)[0]);
  auto iterator(std::next(configuration.begin(), index));
  const auto point(*iterator);
  traits::configuration_manipulation<Configuration>::erase(configuration, iterator);
  return point;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_REMOVE_RANDOM_POINT
