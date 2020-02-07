#ifndef INCLUDE_PPJSDM_RIPPP
#define INCLUDE_PPJSDM_RIPPP

#include <Rcpp.h>
#include <Rinternals.h>

#include "../configuration/configuration_manipulation.hpp"

#include <utility> // std::move

namespace ppjsdm {

template<typename Configuration, typename Window, typename Intensity, typename UpperBound>
inline auto simulate_inhomogeneous_ppp(const Window& window, const Intensity& log_normalised_intensity, const UpperBound& upper_bound, R_xlen_t number_types) {
  const auto volume(window.volume());

  Configuration configuration{};
  // TODO: Might want to .reserve() if Configuration has the member function.
  for(R_xlen_t i(0); i < number_types; ++i) {
    const auto max_points_to_add(R::rpois(volume * static_cast<double>(upper_bound[i])));
    for(R_xlen_t j(0); j < max_points_to_add; ++j) {
      const auto sample(window.sample(i));
      const auto log_normalized_lambda_sample(log_normalised_intensity(sample));
      if(log_normalized_lambda_sample > 0) {
        Rcpp::stop("Did not correctly normalise lambda");
      }
      if(exp_rand() + log_normalized_lambda_sample > 0) {
        add_point(configuration, std::move(sample));
      }
    }
  }
  return configuration;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RIPPP
