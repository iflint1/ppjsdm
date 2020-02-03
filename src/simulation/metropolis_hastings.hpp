#ifndef INCLUDE_PPJSDM_METROPOLIS_HASTINGS
#define INCLUDE_PPJSDM_METROPOLIS_HASTINGS

#include <Rcpp.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "../configuration/configuration_manipulation.hpp"
#include "../utility/size_t.hpp"

namespace ppjsdm {

template<typename Configuration, typename Model, typename Window>
inline auto simulate_metropolis_hastings(const Model& model, const Window& window, unsigned long long int steps, R_xlen_t number_types) {
  constexpr double prob(0.5);
  const auto precomputed_constant((1 - prob) * window.volume() * number_types / prob);

  // TODO: Preallocate with a rough estimate of final size?
  // TODO: Start from non-empty configuration?
  Configuration points{};
  using size_t = size_t<Configuration>;
  size_t points_size(0);

  while(steps-- != 0) {
    if(unif_rand() <= prob) {
      const R_xlen_t random_type(Rcpp::sample(number_types, 1, false, R_NilValue, false)[0]);
      const auto point(window.sample(random_type));

      // TODO: You can avoid taking the exp by reorganising the ratio, and sampling an exponential r.v. instead.
      const auto papangelou(model.compute_papangelou(point, number_types, points));
      const auto birth_ratio(papangelou * precomputed_constant / (1 + points_size));

      // Use C++ short-circuiting
      if(birth_ratio >= 1 || unif_rand() <= birth_ratio) {
        add_point(points, point);
        ++points_size;
      }
    } else if(points_size != 0) {
      const auto saved_point(remove_random_point(points));

      const auto papangelou(model.compute_papangelou(saved_point, number_types, points));
      const auto death_ratio(points_size / (precomputed_constant * papangelou));
      // Use C++ short-circuiting
      if(death_ratio >= 1 || unif_rand() <= death_ratio) {
        --points_size;
      } else {
        add_point(points, saved_point);
      }
    }
  }
  return points;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_METROPOLIS_HASTINGS
