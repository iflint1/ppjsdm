#ifndef INCLUDE_PPJSDM_BIRTH_DEATH
#define INCLUDE_PPJSDM_BIRTH_DEATH

#include <Rcpp.h>

#include "../configuration/configuration_manipulation.hpp"
#include "../configuration/get_number_points.hpp"
#include "../point/point_manipulation.hpp"

#include <limits>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>

namespace ppjsdm {

namespace detail {

// Vectorized rate
template<typename Points, typename Configuration, typename Model>
auto rate(const Points& points, const Configuration& configuration, const Model& model, int nthreads) {
  return model.compute_papangelou_vectorized(points, configuration, nthreads);
}

// Non-vectorized rate
template<typename Point, typename Configuration, typename Model>
auto rate(const Point& point, const Configuration& configuration, const Model& model) {
  return model.compute_papangelou(point, configuration);
}


} // namespace detail

template<typename Configuration, typename Generator, typename Window, typename Quadrature,
         typename BirthModel, typename DeathModel>
inline auto generate_birth_death(Generator& generator,
                                 const Window& window,
                                 const Quadrature& quadrature,
                                 R_xlen_t number_types,
                                 double horizon,
                                 const Configuration& starting_configuration,
                                 const BirthModel& birth_model,
                                 const DeathModel& death_model,
                                 int nthreads) {
  // Random distributions
  std::uniform_int_distribution<int> random_type_distribution(0, number_types - 1);
  std::uniform_real_distribution<double> random_uniform_distribution(0., 1.);
  std::exponential_distribution<double> exponential_distribution(1.);

  // Initialise the configuration
  Configuration configuration(starting_configuration);

  // Initialise a vector of the number of points by type
  std::vector<decltype(get_number_points(configuration))> points_by_type{};
  points_by_type.push_back(get_number_points(configuration, number_types));

  // decltype(detail::rate(quadrature, configuration, birth_model, nthreads)) birth_rates(size(quadrature));
  // Configuration relevant_quadrature(quadrature);
  // std::vector<int> relevant_indices(size(quadrature));
  // for(decltype(relevant_indices.size()) i(0); i < relevant_indices.size(); ++i) {
  //   relevant_indices[i] = i;
  // }

  // Time parameter
  double time(0);

  while(true) {
    // const auto new_birth_rates(detail::rate(relevant_quadrature, configuration, birth_model, nthreads));
    // for(decltype(relevant_indices.size()) i(0); i < relevant_indices.size(); ++i) {
    //   birth_rates[relevant_indices[i]] = new_birth_rates[i];
    // }
    // Compute the birth rates of the quadrature points
    const auto birth_rates(detail::rate(quadrature, configuration, birth_model, nthreads));

    // Constant B approximating \int birth(x, configuration) dx
    const double B(std::accumulate(std::begin(birth_rates), std::end(birth_rates), 0.) /
                   static_cast<double>(birth_rates.size()) * window.volume() * static_cast<double>(number_types));
    // const double B(birth_model.get_upper_bound()[0]); // TODO: Upper bound is per-type, and we are only taking the first type. So only works with one type, for now.

    // Vector containing the death rates of each of the points
    const auto death_rates(detail::rate(configuration, configuration, death_model, nthreads));

    // Compute which of Exp(d_1), ..., Exp(d_N), Exp(B) is smallest
    double running_min(std::numeric_limits<double>::infinity());
    using IntType = decltype(death_rates.size());
    IntType smallest(0);
    for(IntType i(0); i <= death_rates.size(); ++i) {
      double E(exponential_distribution(generator));
      if(i == death_rates.size()) {
        E /= B;
        // E /= birth_rate_max;
      } else {
        E /= death_rates[i];
      }

      if(E < running_min) {
        running_min = E;
        smallest = i;
      }
    }

    // The smallest one gives the time increment
    time += running_min;
    if(time > horizon) {
      break;
    }

    // Initialise the vector containing the number of types to that of the previous step
    points_by_type.push_back(points_by_type.back());

    if(smallest == death_rates.size()) { // Event is a birth
      // Compute (approximate) maximum of the birth rates
      auto birth_rate_max = *max_element(std::begin(birth_rates), std::end(birth_rates));
      if(!configuration.empty()) {
        // Actual maximum is often attained at the configuration points (in particular, when there is attraction)
        // So, also compute the birth rates there and adjust max
        const auto birth_rates_at_configuration(detail::rate(configuration, configuration, birth_model, nthreads));
        birth_rate_max = std::max(birth_rate_max, *max_element(std::begin(birth_rates_at_configuration), std::end(birth_rates_at_configuration)));
      }

      while(true) {
        // Accept sample by computing rejection if rejection sampling condition verified
        const auto random_type(random_type_distribution(generator));
        const auto new_point(window.sample(generator, random_type));
        const auto U(random_uniform_distribution(generator));

        if(U * birth_rate_max < detail::rate(new_point, configuration, birth_model)) { // Accept the sample
          // Update which points of the quadrature are relevant (at which we need to update birth rates)
          // relevant_quadrature = Configuration();
          // relevant_indices = std::vector<int>{};
          // for(decltype(quadrature.size()) i(0); i < quadrature.size(); ++i) {
          //   decltype(birth_model.get_saturation()) close(0);
          //   bool is_relevant(true);
          //   for(decltype(size(configuration)) j(0); j < size(configuration); ++j) {
          //     if((get_x(new_point) - get_x(quadrature[i])) * (get_x(new_point) - get_x(quadrature[i])) +
          //        (get_y(new_point) - get_y(quadrature[i])) * (get_y(new_point) - get_y(quadrature[i])) >
          //        (get_x(quadrature[i]) - get_x(configuration[j])) * (get_x(quadrature[i]) - get_x(configuration[j])) +
          //        (get_y(quadrature[i]) - get_y(configuration[j])) * (get_y(quadrature[i]) - get_y(configuration[j]))) {
          //       ++close;
          //       if(close > birth_model.get_saturation()) {
          //         is_relevant = false;
          //         break;
          //       }
          //     }
          //   }
          //   if(is_relevant) {
          //     relevant_indices.emplace_back(i);
          //     add_point(relevant_quadrature, quadrature[i]);
          //   }
          // }
          add_point(configuration, new_point);
          ++points_by_type.back()[random_type];
          break;
        }
      }
    } else { // Event is a death
      const auto removed_point_iterator(std::begin(configuration) + smallest);
      // relevant_quadrature = Configuration();
      // relevant_indices = std::vector<int>{};
      // for(decltype(quadrature.size()) i(0); i < quadrature.size(); ++i) {
      //   decltype(death_model.get_saturation()) close(0);
      //   bool is_relevant(true);
      //   for(decltype(size(configuration)) j(0); j < size(configuration); ++j) {
      //     if((get_x(*removed_point_iterator) - get_x(quadrature[i])) * (get_x(*removed_point_iterator) - get_x(quadrature[i])) +
      //        (get_y(*removed_point_iterator) - get_y(quadrature[i])) * (get_y(*removed_point_iterator) - get_y(quadrature[i])) >
      //        (get_x(quadrature[i]) - get_x(configuration[j])) * (get_x(quadrature[i]) - get_x(configuration[j])) +
      //        (get_y(quadrature[i]) - get_y(configuration[j])) * (get_y(quadrature[i]) - get_y(configuration[j]))) {
      //       ++close;
      //       if(close >= death_model.get_saturation()) {
      //         is_relevant = false;
      //         break;
      //       }
      //     }
      //   }
      //   if(is_relevant) {
      //     relevant_indices.emplace_back(i);
      //     add_point(relevant_quadrature, quadrature[i]);
      //   }
      // }
      remove_point_by_iterator(configuration, removed_point_iterator);
      --points_by_type.back()[get_type(*removed_point_iterator)];
    }
  }
  return std::make_pair(configuration, points_by_type);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_BIRTH_DEATH
