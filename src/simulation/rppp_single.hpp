#ifndef INCLUDE_PPJSDM_RPPP_SINGLE
#define INCLUDE_PPJSDM_RPPP_SINGLE

#include <Rinternals.h>
#include <Rmath.h>

#include "rbinomialpp_single.hpp"

#include "../utility/window.hpp"

#include <type_traits> // std::remove_cv_t, std::remove_reference_t
#include <vector> // std::vector

namespace ppjsdm {

template<typename Configuration, typename Vector>
inline auto rppp_single(const Window& window, const Vector& lambda) {
  const auto volume(window.volume());
  const auto number_types(lambda.size());
  std::vector<R_xlen_t> number_points(number_types);

  using volume_t = std::remove_cv_t<std::remove_reference_t<decltype(volume)>>;
  for(R_xlen_t type(0); type < number_types; ++type) {
    const auto points_to_add(R::rpois(volume * static_cast<volume_t>(lambda[type])));
    number_points[type] = static_cast<R_xlen_t>(points_to_add);
  }

  return rbinomialpp_single<Configuration>(window, number_points);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RPPP_SINGLE
