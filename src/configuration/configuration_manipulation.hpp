#ifndef INCLUDE_PPJSDM_CONFIGURATION_MANIPULATION
#define INCLUDE_PPJSDM_CONFIGURATION_MANIPULATION

#include <Rcpp.h>
#include <Rinternals.h>

#include <algorithm> // std::find
#include <iterator> // std::iterator_traits, std::next
#include <utility> // std::forward

namespace ppjsdm {
namespace traits {

template<typename Configuration>
struct configuration_manipulation_defaults {
  template<typename Point>
  static inline auto add_point(Configuration& configuration, Point&& point) {
    return configuration.push_back(std::forward<Point>(point));
  }

  template<typename... Args>
  static inline void emplace_point(Configuration& configuration, Args... args) {
    configuration.emplace_back(args...);
  }

  static inline auto size(const Configuration& configuration) {
    return configuration.size();
  }

  static inline bool empty(const Configuration& configuration) {
    return configuration.empty();
  }

  static inline auto remove_point_by_index(Configuration& configuration, R_xlen_t index) {
    auto iterator(std::next(configuration.begin(), index));
    const auto point(*iterator);
    configuration.erase(iterator);
    return point;
  }
};

template<typename Configuration>
struct configuration_manipulation: public configuration_manipulation_defaults<Configuration> {};

} // namespace traits

template<typename Configuration, typename Point>
inline void add_point(Configuration& configuration, Point&& point) {
  traits::configuration_manipulation<Configuration>::add_point(configuration, std::forward<Point>(point));
}

template<typename Configuration, typename... Args>
inline void emplace_point(Configuration& configuration, Args... args) {
  traits::configuration_manipulation<Configuration>::emplace_point(configuration, args...);
}

template<typename Configuration>
inline auto remove_point_by_index(Configuration& configuration, R_xlen_t index) {
  return traits::configuration_manipulation<Configuration>::remove_point_by_index(configuration, index);
}

template<typename Configuration>
inline auto size(const Configuration& configuration) {
  return traits::configuration_manipulation<Configuration>::size(configuration);
}

template<typename Configuration>
inline auto remove_random_point(Configuration& configuration) {
  using difference_type = typename std::iterator_traits<decltype(configuration.begin())>::difference_type;
  const difference_type index(Rcpp::sample(size(configuration), 1, false, R_NilValue, false)[0]);
  return remove_point_by_index(configuration, index);
}

template<typename Configuration>
inline bool empty(const Configuration& configuration) {
  return traits::configuration_manipulation<Configuration>::empty(configuration);
}

// TODO: Does not work well with configuration wrapper class.
template <typename Configuration, typename Point>
inline bool remove_point(Configuration& configuration, const Point& point) {
  auto to_be_erased(std::find(configuration.begin(), configuration.end(), point));
  if(to_be_erased != configuration.end()) {
    configuration.erase(to_be_erased);
    return true;
  }
  return false;
}


} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CONFIGURATION_MANIPULATION
