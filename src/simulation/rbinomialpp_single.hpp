#ifndef INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
#define INCLUDE_PPJSDM_RBINOMIALPP_SINGLE

#include <Rinternals.h>

#include <type_traits> // std::remove_cv_t, std::remove_reference_t

#include "../configuration/configuration_manipulation.hpp"
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
  for(R_xlen_t j(0); j < number_types; ++j) {
    for(filling_t filling(0); filling < n[j]; ++filling) {
      *iterator++ = window.sample(j);
    }
  }

  return configuration;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
