#ifndef INCLUDE_PPJSDM_GET_NUMBER_POINTS
#define INCLUDE_PPJSDM_GET_NUMBER_POINTS

#include "../configuration/configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"

#include <algorithm> // std::max, std::max_element
#include <vector> // std::vector

namespace ppjsdm {

template<typename Configuration>
inline auto get_number_points(const Configuration& configuration,
                              size_t<Configuration> number_types) {
  using size_t = size_t<Configuration>;
  std::vector<size_t> result(number_types);
  for(const auto& point: configuration) {
    const decltype(result.size()) type(get_type(point));
    ++result[type];
  }
  return result;
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
inline auto get_number_points_in_most_numerous_type(const Configuration& configuration) {
  const auto points_by_type(get_number_points(configuration));
  if(!points_by_type.empty()) {
    const auto max(*std::max_element(points_by_type.begin(), points_by_type.end()));
    return max;
  } else {
    return static_cast<typename decltype(points_by_type)::value_type>(0);
  }
}

template<typename Configuration, typename... Configurations>
inline auto get_number_points_in_most_numerous_type(const Configuration& configuration, Configurations&... configurations) {
  return std::max(get_number_points_in_most_numerous_type(configuration), get_number_points_in_most_numerous_type(configurations...));
}

template<typename Configuration>
inline auto get_number_types(const Configuration& configuration) {
  decltype(get_type(*configuration.begin())) max_type(-1);
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
