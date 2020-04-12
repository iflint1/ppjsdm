#ifndef INCLUDE_PPJSDM_RPPP_SINGLE
#define INCLUDE_PPJSDM_RPPP_SINGLE

#include <Rinternals.h>
#include <Rmath.h>

#include "rbinomialpp_single.hpp"

#include <vector> // std::vector

namespace ppjsdm {

template<typename Configuration, typename Lambda>
inline auto rppp_single(const Window& window, const Lambda& lambda, R_xlen_t number_types) {
  const auto volume(window.volume());
  std::vector<R_xlen_t> number_points(number_types);
  R_xlen_t total_number(0);
  for(R_xlen_t j(0); j < number_types; ++j) {
    const auto points_to_add(R::rpois(volume * static_cast<double>(lambda[j])));
    number_points[j] = points_to_add;
    total_number += points_to_add;
  }

  return rbinomialpp_single<Configuration>(window, number_points, number_types, total_number);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RPPP_SINGLE
