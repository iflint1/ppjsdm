#ifndef INCLUDE_PPJSDM_GET_NUMBER_POINTS
#define INCLUDE_PPJSDM_GET_NUMBER_POINTS

#include "configuration_manipulation.h"
#include "../point/point_manipulation.h"

#include <type_traits> // std::remove_const

namespace ppjsdm {

template<typename Configuration>
inline auto get_number_points(const Configuration& configuration, int type) {
  const auto configuration_size(size(configuration));
  using size_t = std::remove_const_t<decltype(configuration_size)>;
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
  return size(configuration);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_GET_NUMBER_POINTS
