#ifndef INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
#define INCLUDE_PPJSDM_RBINOMIALPP_SINGLE

#include <Rinternals.h>

#include "../utility/point_manipulation.h"

namespace ppjsdm {

template<typename Configuration, typename W, typename T>
inline auto rbinomialpp_single(const W& window, const T& n, R_xlen_t point_types, R_xlen_t total_number) {
  Configuration configuration(total_number);

  R_xlen_t fill(0);
  for(R_xlen_t j(0); j < point_types; ++j) {
    const auto points_to_add(n[j]);
    for(R_xlen_t i(0); i < points_to_add; ++i) {
      const auto sample(window.sample());
      configuration[fill + i] = Marked_point(sample.first, sample.second, j);
    }
    fill += points_to_add;
  }

  return configuration;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
