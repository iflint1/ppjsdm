#ifndef INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
#define INCLUDE_PPJSDM_RBINOMIALPP_SINGLE

#include <Rinternals.h>

#include "../configuration/configuration_manipulation.hpp"
#include "../utility/window_utilities.hpp"

namespace ppjsdm {

template<typename Configuration, typename N>
inline auto rbinomialpp_single(const Window& window, const N& n, R_xlen_t number_types, size_t<Configuration> total_number) {
  Configuration configuration(total_number);

  auto iterator(configuration.begin());
  for(R_xlen_t j(0); j < number_types; ++j) {
    // TODO: Why doesn't this work with auto? Type seems to be const-something for some reason.
    R_xlen_t points_to_add(n[j]);
    while(points_to_add-- != 0) {
      *iterator++ = window.sample(j);
    }
  }

  return configuration;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
