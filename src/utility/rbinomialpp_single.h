#ifndef INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
#define INCLUDE_PPJSDM_RBINOMIALPP_SINGLE

#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration_wrapper.h"
#include "point_manipulation.h"

namespace ppjsdm {

template<typename Configuration, typename W, typename T>
inline auto rbinomialpp_single(const W& window, const T& n, R_xlen_t point_types, R_xlen_t total_number) {
  Configuration r_configuration(total_number);

  R_xlen_t fill(0);
  for(R_xlen_t j(0); j < point_types; ++j) {
    const auto points_to_add(n[j]);
    for(R_xlen_t i(0); i < points_to_add; ++i) {
      const auto sample(window.sample());
      r_configuration[fill + i] = Marked_point(sample.first, sample.second, j);
    }
    fill += points_to_add;
  }

  return r_configuration;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
