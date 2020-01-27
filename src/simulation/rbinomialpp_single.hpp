#ifndef INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
#define INCLUDE_PPJSDM_RBINOMIALPP_SINGLE

#include <Rinternals.h>

#include "../utility/size_t.hpp"

namespace ppjsdm {

template<typename Configuration, typename Window, typename N>
inline auto rbinomialpp_single(const Window& window, const N& n, R_xlen_t number_types, size_t<Configuration> total_number) {
  Configuration configuration(total_number);
  using size_t = size_t<Configuration>;

  size_t fill(0);
  for(R_xlen_t j(0); j < number_types; ++j) {
    const auto points_to_add = static_cast<size_t>(n[j]);
    for(size_t i(0); i < points_to_add; ++i) {
      configuration[fill + i] = window.sample(j);
    }
    fill += points_to_add;
  }

  return configuration;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RBINOMIALPP_SINGLE