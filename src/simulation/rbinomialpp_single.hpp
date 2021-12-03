#ifndef INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
#define INCLUDE_PPJSDM_RBINOMIALPP_SINGLE

#include <Rinternals.h>

#include <type_traits> // std::remove_cv_t, std::remove_reference_t

#include "../utility/window.hpp"

namespace ppjsdm {

template<typename Configuration, typename Vector>
inline auto rbinomialpp_single(const Window& window,
                               const Vector& n,
                               R_xlen_t number_types,
                               size_t<Configuration> total_number) {
  Configuration configuration(total_number);

  auto iterator(configuration.begin());
  using filling_t = std::remove_cv_t<std::remove_reference_t<decltype(n[0])>>;
  for(R_xlen_t type(0); type < number_types; ++type) {
    for(filling_t filling(0); filling < n[type]; ++filling) {
      *iterator++ = window.sample(type);
    }
  }

  return configuration;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
