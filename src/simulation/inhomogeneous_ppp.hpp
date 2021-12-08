#ifndef INCLUDE_PPJSDM_RIPPP
#define INCLUDE_PPJSDM_RIPPP

#include <Rcpp.h>
#include <Rinternals.h>

#include "../configuration/configuration_manipulation.hpp"
#include "../utility/window.hpp"

#include <utility> // std::move

namespace ppjsdm {

template<typename Configuration, typename Intensity, typename UpperBound>
inline auto simulate_inhomogeneous_ppp(const Window& window,
                                       const Intensity& log_normalised_intensity,
                                       const UpperBound& upper_bound) {
  const auto volume(window.volume());
  const auto number_types(upper_bound.size());

  Configuration configuration{};

  using volume_t = std::remove_cv_t<std::remove_reference_t<decltype(volume)>>;
  for(R_xlen_t type(0); type < number_types; ++type) {
    const R_xlen_t max_points_to_add(R::rpois(volume * static_cast<volume_t>(upper_bound[type])));
    reserve_if_possible(configuration, configuration.size() + max_points_to_add);
    for(R_xlen_t j(0); j < max_points_to_add; ++j) {
      const auto sample(window.sample(type));
      const auto log_normalized_lambda_sample(log_normalised_intensity(sample));
      if(log_normalized_lambda_sample > 0) {
        Rcpp::stop("Did not correctly normalise the intensity.");
      }
      // Use C++ short-circuiting
      if(log_normalized_lambda_sample >= 0 || exp_rand() + log_normalized_lambda_sample >= 0) {
        add_point(configuration, std::move(sample));
      }
    }
  }
  return configuration;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RIPPP
