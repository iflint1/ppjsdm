#ifndef INCLUDE_PPJSDM_CONFIGURATION_MANIPULATION
#define INCLUDE_PPJSDM_CONFIGURATION_MANIPULATION

#include "configuration_wrappers.h"

namespace ppjsdm {
namespace traits {

template<typename Configuration>
struct configuration_manipulation_defaults {
  template<typename... Args>
  static inline void emplace_point(Configuration& configuration, Args... args) {
    configuration.emplace_back(args...);
  }

  static inline auto size(const Configuration& configuration) {
    return configuration.size();
  }
};

template<typename Configuration>
struct configuration_manipulation: public configuration_manipulation_defaults<Configuration> {};

template<>
struct configuration_manipulation<Configuration_wrapper>:
public configuration_manipulation_defaults<Configuration_wrapper> {
  using Configuration = Configuration_wrapper;
  template<typename... Args>
  static inline void emplace_point(Configuration& configuration, Args... args) {
    configuration.push_back(args...);
  }
};

} // namespace traits

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
