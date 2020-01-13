#ifndef INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
#define INCLUDE_PPJSDM_RBINOMIALPP_SINGLE

#include <Rinternals.h>

#include "../point/point_manipulation.h"

namespace ppjsdm {

template<typename Configuration, typename Window, typename N>
inline auto rbinomialpp_single(const Window& window, const N& n, R_xlen_t point_types, R_xlen_t total_number) {
  Configuration configuration(total_number);

  R_xlen_t fill(0);
  for(R_xlen_t j(0); j < point_types; ++j) {
    const auto points_to_add(n[j]);
    for(R_xlen_t i(0); i < points_to_add; ++i) {
      configuration[fill + i] = window.sample(j);
    }
    fill += points_to_add;
  }

  return configuration;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
