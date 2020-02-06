#ifndef INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST
#define INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST

#include <Rcpp.h>
#include <Rinternals.h>

#include "../configuration/configuration_manipulation.hpp"
#include "../simulation/inhomogeneous_ppp.hpp"
#include "../utility/backwards_markov_chain.hpp"

#include <tuple> // std::pair

namespace ppjsdm {
namespace detail {

// This is an efficient implementation of the generic CFTP algorithm from Moller and Waagepetersen's book, cf. p. 230 therein.
// The point process is neither assumed to be attractive nor repulsive, see also ``Perfect simulation of spatial
// point processes using dominated coupling from the past with application to a multiscale area-interaction point process''
// by Ambler and Silverman for a similar (but less general) technique.
template<typename Chain, typename Model, typename Window>
inline auto update_LU_and_check_coalescence(Chain& chain, const Model& model, const Window& window) {
  auto points_not_in_L(chain.get_last_configuration()); // U starts with the end value of Z
  using Configuration = decltype(points_not_in_L);
  Configuration L{}; // L is an empty point process
  chain.iterate_forward_in_time([&points_not_in_L, &L, &model, &window](auto&& point, auto uniform_mark) {
    // TODO: Might be able to reorganise, take log() and get speed-up.
    const auto alpha_max(model.compute_log_alpha_max(window, point, L, points_not_in_L));
    if(alpha_max > 1) {
      Rcpp::stop("Did not correctly normalize the Papangelou intensity");
    }
    if(alpha_max > uniform_mark) {
      const auto alpha_min(model.compute_alpha_min(window, point, L, points_not_in_L));
      if(alpha_min > 1) {
        Rcpp::stop("Did not correctly normalize the Papangelou intensity");
      }
      add_point(alpha_min > uniform_mark ? L : points_not_in_L, std::forward<decltype(point)>(point));
    }
  }, [&points_not_in_L, &L](auto&& point) {
    if(!remove_point(points_not_in_L, point)) {
      remove_point(L, std::forward<decltype(point)>(point));
    }
  });
  return std::pair<bool, Configuration>(empty(points_not_in_L), L);
}

} // namespace detail

template<typename Configuration, typename Model, typename Window>
inline auto simulate_coupling_from_the_past(const Model& model, const Window& window, R_xlen_t number_types) {
  const auto normalised_dominating_intensity(model.get_normalised_dominating_intensity(window));
  const auto intensity_upper_bound(normalised_dominating_intensity.get_upper_bound());
  Backwards_Markov_chain<Configuration> Z(simulate_inhomogeneous_ppp<Configuration>(window, normalised_dominating_intensity, intensity_upper_bound, number_types));

  const auto integral_of_dominating_intensity(normalised_dominating_intensity.get_integral());
  const auto T0(Z.extend_until_T0(integral_of_dominating_intensity, window, number_types));
  while(true) {
    const auto coalescence(detail::update_LU_and_check_coalescence(Z, model, window));
    if(coalescence.first) {
      return coalescence.second;
    }
    R_CheckUserInterrupt();
    Z.extend_backwards(T0, integral_of_dominating_intensity, window, number_types);
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST
