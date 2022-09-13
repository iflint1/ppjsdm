#ifndef INCLUDE_PPJSDM_RIPPP
#define INCLUDE_PPJSDM_RIPPP

#include <Rcpp.h>
#include <Rinternals.h>

#include "../configuration/configuration_manipulation.hpp"
#include "../utility/window.hpp"

#include <random> // Distribution
#include <utility> // std::move

namespace ppjsdm {

template<typename Configuration, typename Generator, typename Intensity, typename UpperBound>
inline auto simulate_inhomogeneous_ppp(Generator& generator,
                                       const Window& window,
                                       const Intensity& log_normalised_intensity,
                                       const UpperBound& upper_bound) {
  const auto volume(window.volume());
  const auto number_types(upper_bound.size());

  Configuration configuration{};

  // Probability distribution
  std::exponential_distribution<double> exponential_distribution(1);

  using volume_t = std::remove_cv_t<std::remove_reference_t<decltype(volume)>>;
  using number_types_t = std::remove_cv_t<std::remove_reference_t<decltype(number_types)>>;
  for(number_types_t type(0); type < number_types; ++type) {
    const R_xlen_t max_points_to_add(R::rpois(volume * static_cast<volume_t>(upper_bound[type])));
    reserve_if_possible(configuration, configuration.size() + max_points_to_add);
    for(R_xlen_t j(0); j < max_points_to_add; ++j) {
      const auto sample(window.sample(generator, type));
      const auto log_normalized_lambda_sample(log_normalised_intensity(sample));
      if(log_normalized_lambda_sample > 0) {
        Rcpp::stop("Did not correctly normalise the intensity.");
      }
      // Use C++ short-circuiting
      if(log_normalized_lambda_sample >= 0 || exponential_distribution(generator) + log_normalized_lambda_sample >= 0) {
        add_point(configuration, std::move(sample));
      }
    }
  }
  return configuration;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RIPPP
