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
// The point process is neither assumed to be attractive nor repulsive, see also p. 361 of ``A primer on perfect simulation
// for spatial point processes'' by Berthelsen and Moller, and ``Perfect simulation of spatial
// point processes using dominated coupling from the past with application to a multiscale area-interaction point process''
// by Ambler and Silverman for a similar (but less general) idea.
template<typename Chain, typename Model, typename Window>
inline auto update_LU_and_check_coalescence(Chain& chain, const Model& model, const Window& window) {
  auto points_not_in_L(chain.get_last_configuration()); // U starts with the end value of Z
  using Configuration = decltype(points_not_in_L);
  Configuration L{}; // L is an empty point process
  chain.iterate_forward_in_time([&points_not_in_L, &L, &model, &window](auto&& point, auto exp_mark) {
    const auto log_alpha_max(model.compute_log_alpha_max(window, point, L, points_not_in_L));
    if(log_alpha_max > 0) {
      Rcpp::stop("Did not correctly normalize the Papangelou intensity");
    }
    if(log_alpha_max + exp_mark > 0) {
      const auto log_alpha_min(model.compute_log_alpha_min(window, point, L, points_not_in_L));
      if(log_alpha_min > 0) {
        Rcpp::stop("Did not correctly normalize the Papangelou intensity");
      }
      add_point(log_alpha_min + exp_mark > 0 ? L : points_not_in_L, std::forward<decltype(point)>(point));
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
  auto Z(make_backwards_markov_chain<Configuration>(model.get_normalised_dominating_intensity(window), number_types));

  Z.extend_until_T0();
  while(true) {
    const auto coalescence(detail::update_LU_and_check_coalescence(Z, model, window));
    if(coalescence.first) {
      return coalescence.second;
    }
    R_CheckUserInterrupt();
    Z.extend_backwards(Z.size());
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST
