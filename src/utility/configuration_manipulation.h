#ifndef INCLUDE_PPJSDM_CONFIGURATION_MANIPULATION
#define INCLUDE_PPJSDM_CONFIGURATION_MANIPULATION

#include "configuration_wrappers.h"

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

  template<typename Iterator>
  static inline auto erase(Configuration& configuration, Iterator it) {
    return configuration.erase(it);
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

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CONFIGURATION_MANIPULATION
