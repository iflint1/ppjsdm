#ifndef INCLUDE_PPJSDM_METROPOLIS_HASTINGS
#define INCLUDE_PPJSDM_METROPOLIS_HASTINGS

#include <Rcpp.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "../configuration/configuration_manipulation.hpp"

namespace ppjsdm {

template<typename Configuration, typename Model, typename Window>
inline auto simulate_metropolis_hastings(const Model& model, const Window& window, unsigned long long int steps, R_xlen_t number_types) {
  constexpr double prob(0.5);
  const double precomputed_constant((1 - prob) * window.volume() * static_cast<double>(number_types) / prob);

  // TODO: Preallocate with a rough estimate of final size?
  // TODO: Start from non-empty configuration?
  Configuration points{};

  while(steps-- != 0) {
    if(unif_rand() <= prob) {
      const R_xlen_t random_type(Rcpp::sample(number_types, 1, false, R_NilValue, false)[0]);
      const auto point(window.sample(random_type));

      // TODO: You can avoid taking the exp by reorganising the ratio, and sampling an exponential r.v. instead.
      const double papangelou(model.compute_papangelou(point, points));
      const double birth_ratio(papangelou * precomputed_constant / (1 + size(points)));

      // Use C++ short-circuiting
      if(birth_ratio >= 1. || unif_rand() <= birth_ratio) {
        add_point(points, point);
      }
    } else if(!empty(points)) {
      const auto saved_point(remove_random_point(points));

      const double papangelou(model.compute_papangelou(saved_point, points));
      const double death_ratio(static_cast<double>(size(points)) / (precomputed_constant * papangelou));
      // Use C++ short-circuiting
      if(death_ratio < 1. && unif_rand() > death_ratio) {
        add_point(points, saved_point);
      }
    }
  }
  return points;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_METROPOLIS_HASTINGS
