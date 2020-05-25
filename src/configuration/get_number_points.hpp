#ifndef INCLUDE_PPJSDM_GET_NUMBER_POINTS
#define INCLUDE_PPJSDM_GET_NUMBER_POINTS

#include "../configuration/configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"

#include <vector> // std::vector

namespace ppjsdm {

template<typename Configuration>
inline auto get_number_points(const Configuration& configuration, int type) {
  using size_t = size_t<Configuration>;
  size_t total(0);
  for(const auto& point: configuration) {
    if(get_type(point) == type) {
      ++total;
    }
  }
  return total;
}

template<typename Configuration>
inline auto get_number_points(const Configuration& configuration) {
  using size_t = size_t<Configuration>;
  std::vector<size_t> result{};
  for(const auto& point: configuration) {
    const decltype(result.size()) type(get_type(point));
    while(type >= result.size()) {
      result.emplace_back(0);
    }
    ++result[type];
  }
  return result;
}

template<typename Configuration>
inline auto get_number_types(const Configuration& configuration) {
  int max_type(-1);
  for(const auto& point: configuration) {
    const auto type(get_type(point));
    if(type > max_type) {
      max_type = type;
    }
  }
  return max_type + 1;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_GET_NUMBER_POINTS
