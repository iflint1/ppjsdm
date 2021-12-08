#ifndef INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
#define INCLUDE_PPJSDM_RBINOMIALPP_SINGLE

#include <Rinternals.h>

#include <type_traits> // std::remove_cv_t, std::remove_reference_t

#include "../utility/sum.hpp"
#include "../utility/window.hpp"

namespace ppjsdm {

template<typename Configuration, typename Vector>
inline auto rbinomialpp_single(const Window& window,
                               const Vector& n) {
  const auto total_number(sum<R_xlen_t>(n));
  Configuration configuration(total_number);

  auto iterator(configuration.begin());

  const auto number_types(n.size());

  using filling_t = std::remove_cv_t<std::remove_reference_t<decltype(n[0])>>;
  using type_t = decltype(n.size());
  for(type_t type(0); type < number_types; ++type) {
    for(filling_t filling(0); filling < n[type]; ++filling) {
      *iterator++ = window.sample(type);
    }
  }

  return configuration;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
