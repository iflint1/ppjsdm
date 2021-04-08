#ifndef INCLUDE_PPJSDM_METROPOLIS_HASTINGS
#define INCLUDE_PPJSDM_METROPOLIS_HASTINGS

#include <Rcpp.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "../configuration/configuration_manipulation.hpp"
#include "../utility/approximate_draw.hpp"

namespace ppjsdm {

template<typename Configuration, typename Model>
inline auto simulate_metropolis_hastings(const Model& model, const Window& window, unsigned long long int steps, R_xlen_t number_types) {
  // Start from a rough approximate draw and go from there.
  const auto points(approximate_draw<Configuration>(model));

  return simulate_metropolis_hastings(model, window, steps, number_types, points);
}

template<typename Configuration, typename Model>
inline auto simulate_metropolis_hastings(const Model& model, const Window& window, unsigned long long int steps, R_xlen_t number_types, const Configuration& starting_configuration) {
  constexpr double birth_probability(0.5);
  const double precomputed_constant((1 - birth_probability) * window.volume() * static_cast<double>(number_types) / birth_probability);

  Configuration points(starting_configuration);

  for(decltype(steps) k(0); k < steps; ++k) {
    if(unif_rand() <= birth_probability) {
      const R_xlen_t random_type(Rcpp::sample(number_types, 1, false, R_NilValue, false)[0]);
      const auto point(window.sample(random_type));

      // TODO: You can avoid taking the exp by reorganising the ratio, and sampling an exponential r.v. instead.
      const double papangelou(model.compute_papangelou(point, points));
      const double birth_ratio(papangelou * precomputed_constant / (1. + static_cast<double>(size(points))));

      // Use C++ short-circuiting
      if(birth_ratio >= 1. || unif_rand() <= birth_ratio) {
        add_point(points, point);
      }
    } else if(!empty(points)) {
      const double current_size(size(points));
      const auto point(remove_random_point(points));

      const double papangelou(model.compute_papangelou(point, points));
      double death_ratio(current_size / (precomputed_constant * papangelou));

      // Use C++ short-circuiting
      if(death_ratio < 1. && unif_rand() > death_ratio) {
        add_point(points, point);
      }
    }
  }
  return points;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_METROPOLIS_HASTINGS
