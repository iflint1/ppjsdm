#ifndef INCLUDE_PPJSDM_METROPOLIS_HASTINGS
#define INCLUDE_PPJSDM_METROPOLIS_HASTINGS

#include <Rcpp.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "../configuration/configuration_manipulation.hpp"
#include "../utility/approximate_draw.hpp"
#include "../utility/timer.hpp"

#include <cmath> // std::floor
#include <iterator> // std::next
#include <vector> // std::vector

namespace ppjsdm {

template<typename Configuration, typename Model>
inline auto simulate_metropolis_hastings(const Model& model, unsigned long long int steps) {
  // Start from a rough approximate draw and go from there.
  const auto starting_configuration(approximate_draw<Configuration>(model));

  return simulate_metropolis_hastings(model, steps, starting_configuration);
}

template<typename Configuration, typename Model>
inline auto simulate_metropolis_hastings(const Model& model, unsigned long long int steps, const Configuration& starting_configuration) {
  constexpr double birth_probability(0.5);
  const auto number_types(model.get_number_types());
  const auto window(model.get_window());
  const double precomputed_constant((1 - birth_probability) / birth_probability * window.volume() * static_cast<double>(number_types));

  Configuration points(starting_configuration);

  // Set up timer to regularly check for user interruption
  PreciseTimer timer{};
  timer.set_current();

  for(decltype(steps) k(0); k < steps; ++k) {
    // Allow user interruption at regular intervals
    const auto t(timer.get_total_time().count());
    if((static_cast<int>(std::floor(t * 10)) % 10) == 0) {
      Rcpp::checkUserInterrupt();
    }
    if(unif_rand() <= birth_probability) {
      const R_xlen_t random_type(Rcpp::sample(number_types, 1, false, R_NilValue, false)[0]);
      const auto point(window.sample(random_type));

      const double papangelou(model.compute_papangelou(point, points));
      const double birth_ratio(papangelou * precomputed_constant / (1. + static_cast<double>(size(points))));

      // Use C++ short-circuiting
      if(birth_ratio >= 1. || birth_ratio >= unif_rand()) {
        add_point(points, point);
      }
    } else if(!empty(points)) {
      const double current_size(size(points));
      const auto point(remove_random_point(points));

      const double papangelou(model.compute_papangelou(point, points));
      double death_ratio(current_size / (precomputed_constant * papangelou));

      // Use C++ short-circuiting
      if(death_ratio <= 1. && death_ratio <= unif_rand()) {
        add_point(points, point);
      }
    }
  }
  // for(decltype(steps) k(0); k < steps; ++k) {
  //   if(unif_rand() <= birth_probability) {
  //     const R_xlen_t random_type(Rcpp::sample(number_types, 1, false, R_NilValue, false)[0]);
  //     const auto point(window.sample(random_type));
  //
  //     const double log_papangelou(model.compute_log_papangelou(point, points));
  //     const double log_birth_ratio(log_papangelou + precomputed_constant - std::log(1. + static_cast<double>(size(points))));
  //
  //     // Use C++ short-circuiting
  //     if(log_birth_ratio >= 0. || log_birth_ratio + exp_rand() >= 0) {
  //       add_point(points, point);
  //     }
  //   } else if(!empty(points)) {
  //     const double current_size(size(points));
  //     const auto point(remove_random_point(points));
  //
  //     const double log_papangelou(model.compute_log_papangelou(point, points));
  //     const double log_death_ratio(log_papangelou + precomputed_constant - std::log(current_size));
  //
  //     // Use C++ short-circuiting
  //     if(log_death_ratio >= 0. && log_death_ratio >= exp_rand()) {
  //       add_point(points, point);
  //     }
  //   }
  // }
  return points;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_METROPOLIS_HASTINGS
