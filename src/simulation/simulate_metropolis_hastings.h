#ifndef INCLUDE_PPJSDM_RMULTIGIBBS_SINGLE
#define INCLUDE_PPJSDM_RMULTIGIBBS_SINGLE

#include <Rinternals.h>
#include <Rmath.h>

#include "../utility/configuration_manipulation.h"

namespace ppjsdm {

template<typename Configuration, typename Model, typename Window>
inline auto simulate_metropolis_hastings(const Model& model, const Window& window, R_xlen_t steps, R_xlen_t point_types) {
  const auto volume(window.volume());
  constexpr double prob(0.5);

  // TODO: Preallocate with a rough estimate of final size?
  // TODO: Start from non-empty configuration?
  Configuration points{};
  R_xlen_t total_number(0);

  for(R_xlen_t step(0); step < steps; ++step) {
    const auto u(unif_rand());
    const auto v(unif_rand());
    if(u <= prob) {
      const R_xlen_t type(Rcpp::sample(point_types, 1, false, R_NilValue, false)[0]);
      const auto point(window.sample(type));

      // TODO: You can avoid taking the exp by reorganising the ratio, and sampling an exponential r.v. instead.
      const auto papangelou(model.compute_papangelou(points, point, point_types));
      const auto birth_ratio(papangelou * (1 - prob) * volume * point_types / (prob * (1 + total_number)));

      if(v <= birth_ratio) {
        ppjsdm::add_point(points, point);
        ++total_number;
      }
    } else {
      if(total_number != 0) {
        const auto saved_point(ppjsdm::remove_random_point(points));

        const auto papangelou(model.compute_papangelou(points, saved_point, point_types));
        const auto death_ratio(prob * total_number / (point_types * (1 - prob) * volume * papangelou));
        if(v <= death_ratio) {
          --total_number;
        } else {
          ppjsdm::add_point(points, saved_point);
        }
      }
    }
  }
  return points;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RMULTIGIBBS_SINGLE
