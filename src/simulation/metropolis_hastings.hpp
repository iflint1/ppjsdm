#ifndef INCLUDE_PPJSDM_METROPOLIS_HASTINGS
#define INCLUDE_PPJSDM_METROPOLIS_HASTINGS

#include <Rcpp.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "../configuration/configuration_manipulation.hpp"
#include "../utility/size_t.hpp"

namespace ppjsdm {

template<typename Configuration, typename Model, typename Window>
inline auto simulate_metropolis_hastings(const Model& model, const Window& window, unsigned long long int steps, R_xlen_t point_types) {
  const auto volume(window.volume());
  constexpr double prob(0.5);
  using size_t = size_t<Configuration>;

  // TODO: Preallocate with a rough estimate of final size?
  // TODO: Start from non-empty configuration?
  Configuration points{};
  size_t total_number(0);

  while(steps-- != 0) {
    const auto u(unif_rand());
    const auto v(unif_rand());
    if(u <= prob) {
      const R_xlen_t type(Rcpp::sample(point_types, 1, false, R_NilValue, false)[0]);
      const auto point(window.sample(type));

      // TODO: You can avoid taking the exp by reorganising the ratio, and sampling an exponential r.v. instead.
      const auto papangelou(model.compute_papangelou(points, point, point_types));
      const auto birth_ratio(papangelou * (1 - prob) * volume * point_types / (prob * (1 + total_number)));

      if(v <= birth_ratio) {
        add_point(points, point);
        ++total_number;
      }
    } else {
      if(total_number != 0) {
        const auto saved_point(remove_random_point(points));

        const auto papangelou(model.compute_papangelou(points, saved_point, point_types));
        const auto death_ratio(prob * total_number / (point_types * (1 - prob) * volume * papangelou));
        if(v <= death_ratio) {
          --total_number;
        } else {
          add_point(points, saved_point);
        }
      }
    }
  }
  return points;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_METROPOLIS_HASTINGS
