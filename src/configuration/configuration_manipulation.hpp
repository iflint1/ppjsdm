#ifndef INCLUDE_PPJSDM_CONFIGURATION_MANIPULATION
#define INCLUDE_PPJSDM_CONFIGURATION_MANIPULATION

#include <Rcpp.h>
#include <Rinternals.h>

#include <algorithm> // std::find
#include <iterator> // std::iterator_traits, std::next
#include <utility> // std::forward, std::declval

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

  static inline auto remove_point_by_iterator(Configuration& configuration, typename Configuration::iterator iterator) {
    const auto point(*iterator);
    configuration.erase(iterator);
    return point;
  }

  static inline void reserve_if_possible(Configuration& configuration, typename Configuration::size_type n) {
    configuration.reserve(n);
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
inline auto size(const Configuration& configuration) {
  return traits::configuration_manipulation<Configuration>::size(configuration);
}

template<typename Configuration, typename... Configurations, std::enable_if_t<sizeof...(Configurations) != 0>* = nullptr>
inline auto size(const Configuration& configuration, Configurations&... configurations) {
  return traits::configuration_manipulation<Configuration>::size(configuration) + size(configurations...);
}

template<typename T>
using size_t = decltype(size(std::declval<T>()));

template<typename Configuration>
inline auto remove_point_by_iterator(Configuration& configuration, typename Configuration::iterator iterator) {
  return traits::configuration_manipulation<Configuration>::remove_point_by_iterator(configuration, iterator);
}

template<typename Configuration>
inline auto remove_random_point(Configuration& configuration) {
  using difference_type = typename std::iterator_traits<decltype(configuration.begin())>::difference_type;
  const difference_type index(Rcpp::sample(size(configuration), 1, false, R_NilValue, false)[0]);
  return remove_point_by_iterator(configuration, std::next(configuration.begin(), index));
}

template<typename Configuration>
inline auto random_point(Configuration& configuration) {
  using difference_type = typename std::iterator_traits<decltype(configuration.begin())>::difference_type;
  const difference_type index(Rcpp::sample(size(configuration), 1, false, R_NilValue, false)[0]);
  return std::next(configuration.begin(), index);
}

// TODO: Combine with function above (through a remove_by_index?)
template<typename Generator, typename Configuration>
inline auto remove_random_point(Generator& generator, Configuration& configuration) {
  using difference_type = typename std::iterator_traits<decltype(configuration.begin())>::difference_type;
  std::uniform_int_distribution<difference_type> random_point_distribution(0, size(configuration) - 1);
  const auto index(random_point_distribution(generator));
  return remove_point_by_iterator(configuration, std::next(configuration.begin(), index));
}

template<typename Configuration>
inline bool empty(const Configuration& configuration) {
  return traits::configuration_manipulation<Configuration>::empty(configuration);
}

template <typename Configuration, typename Point>
inline bool remove_point(Configuration& configuration, const Point& point) {
  const auto iterator(std::find(configuration.begin(), configuration.end(), point));
  if(iterator != configuration.end()) {
    remove_point_by_iterator(configuration, iterator);
    return true;
  }
  return false;
}

template <typename Configuration>
inline void reserve_if_possible(Configuration& configuration, typename Configuration::size_type n) {
  traits::configuration_manipulation<Configuration>::reserve_if_possible(configuration, n);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CONFIGURATION_MANIPULATION
