#ifndef INCLUDE_PPJSDM_GET_NUMBER_POINTS
#define INCLUDE_PPJSDM_GET_NUMBER_POINTS

#include "../configuration/configuration_manipulation.h"
#include "../point/point_manipulation.h"
#include "../utility/size_t.h"

#include <vector> // std::vector

namespace ppjsdm {

template<typename Configuration>
inline auto get_number_points(const Configuration& configuration, int type) {
  const auto configuration_size(size(configuration));
  using size_t = size_t<Configuration>;
  size_t total(0);
  for(size_t i(0); i < configuration_size; ++i) {
    if(get_type(configuration[i]) == type) {
      ++total;
    }
  }
  return total;
}

template<typename Configuration>
inline auto get_number_points(const Configuration& configuration) {
  const auto configuration_size(size(configuration));
  using size_t = size_t<Configuration>;
  std::vector<size_t> result{};
  for(size_t i(0); i < configuration_size; ++i) {
    const decltype(result.size()) type(get_type(configuration[i]));
    while(type >= result.size()) {
      result.emplace_back(0);
    }
    ++result[type];
  }
  return result;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_GET_NUMBER_POINTS
