#ifndef INCLUDE_PPJSDM_METROPOLIS_HASTINGS
#define INCLUDE_PPJSDM_METROPOLIS_HASTINGS

#include <Rcpp.h>
#include <RcppThread.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "../configuration/configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"
#include "../utility/approximate_draw.hpp"
#include "../utility/timer.hpp"

#include <cmath> // std::floor
#include <iterator> // std::next
#include <random> // Generator, distribution, etc.
#include <vector> // std::vector

namespace ppjsdm {

template<typename Configuration, typename Generator, typename Model, typename Vector>
inline auto simulate_metropolis_hastings(Generator& generator,
                                         const Model& model,
                                         unsigned long long int steps,
                                         const Vector& only_simulate_these_types,
                                         const Configuration& conditional_configuration) {
  // Start from a rough approximate draw and go from there.
  const auto starting_configuration(approximate_draw<Configuration>(generator, model));

  // Remove points which are not of the right type
  Configuration thinned_starting_configuration{};
  reserve_if_possible(thinned_starting_configuration, size(starting_configuration));
  for(const auto& point: starting_configuration) {
    for(decltype(only_simulate_these_types.size()) j(0); j < only_simulate_these_types.size(); ++j) {
      if(get_type(point) == only_simulate_these_types[j]) {
        add_point(thinned_starting_configuration, point);
        break;
      }
    }
  }

  return simulate_metropolis_hastings(generator,
                                      model,
                                      steps,
                                      thinned_starting_configuration,
                                      only_simulate_these_types,
                                      conditional_configuration);
}

template<typename Configuration, typename Generator, typename Model, typename Vector>
inline auto simulate_metropolis_hastings(Generator& generator,
                                         const Model& model,
                                         unsigned long long int steps,
                                         const Configuration& starting_configuration,
                                         const Vector& only_simulate_these_types,
                                         const Configuration& conditional_configuration) {
  // Probability of drawing tentative births
  constexpr double birth_probability(0.5);

  // Number of types we are simulating from
  const auto number_types(only_simulate_these_types.size());

  // Random distributions
  std::uniform_int_distribution<decltype(only_simulate_these_types.size())> random_type_distribution(0, number_types - 1);
  std::uniform_real_distribution<double> random_uniform_distribution(0, 1);

  // Simulation window
  const auto window(model.get_window());

  // Can precompute this
  const double precomputed_constant((1 - birth_probability) / birth_probability * window.volume() * static_cast<double>(number_types));

  Configuration points(starting_configuration);

  // Set up timer to regularly check for user interruption
  PreciseTimer timer{};
  timer.set_current();

  for(decltype(steps) k(0); k < steps; ++k) {
    // Allow user interruption at regular intervals
    const auto t(timer.get_total_time().count());
    if((static_cast<int>(std::floor(t * 10)) % 10) == 0) {
      RcppThread::checkUserInterrupt();
    }

    if(random_uniform_distribution(generator) <= birth_probability) { // Births
      const auto random_type(only_simulate_these_types[random_type_distribution(generator)]);
      const auto point(window.sample(generator, random_type));

      const double papangelou(model.compute_papangelou(point, points, conditional_configuration));
      const double birth_ratio(papangelou * precomputed_constant / (1. + static_cast<double>(size(points))));

      // Use C++ short-circuiting
      if(birth_ratio >= 1. || birth_ratio >= random_uniform_distribution(generator)) {
        add_point(points, point);
      }
    } else if(!empty(points)) { // Deaths
      const auto current_size = static_cast<double>(size(points));
      const auto point(remove_random_point(generator, points));

      const double papangelou(model.compute_papangelou(point, points, conditional_configuration));
      const double death_ratio(current_size / (precomputed_constant * papangelou));

      // Use C++ short-circuiting
      if(death_ratio <= 1. && death_ratio <= random_uniform_distribution(generator)) {
        add_point(points, point);
      }
    }
  }
  return points;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_METROPOLIS_HASTINGS
