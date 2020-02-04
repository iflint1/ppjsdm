#ifndef INCLUDE_PPJSDM_RIPPP
#define INCLUDE_PPJSDM_RIPPP

#include <Rcpp.h>
#include <Rinternals.h>

#include "../configuration/configuration_manipulation.hpp"

namespace ppjsdm {

template<typename Configuration, typename Window, typename Lambda, typename UpperBound>
inline auto simulate_inhomogeneous_ppp(const Window& window, const Lambda& normalized_lambda, const UpperBound& upper_bound, R_xlen_t number_types) {
  const auto volume(window.volume());

  Configuration configuration{};
  // TODO: Might want to .reserve() if Configuration has the member function.
  for(R_xlen_t i(0); i < number_types; ++i) {
    const auto max_points_to_add(R::rpois(volume * static_cast<double>(upper_bound[i])));
    for(R_xlen_t j(0); j < max_points_to_add; ++j) {
      const auto sample(window.sample(i));
      const auto normalized_lambda_sample(normalized_lambda(sample));
      if(normalized_lambda_sample < 0 || normalized_lambda_sample > 1) {
        Rcpp::stop("Did not correctly normalise lambda");
      }
      // TODO: As usual, might be able to take log on both sides...
      if(unif_rand() < normalized_lambda_sample) {
        add_point(configuration, std::move(sample));
      }
    }
  }
  return configuration;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RIPPP
