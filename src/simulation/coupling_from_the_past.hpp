#ifndef INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST
#define INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST

#include <Rinternals.h>

#include "../simulation/rppp_single.hpp"
#include "../utility/backwards_markov_chain.hpp"

#include <tuple> // std::pair

namespace ppjsdm {
namespace detail {

template<typename Lambda, typename Window>
inline auto get_integral_of_dominating_intensity(const Lambda& dominating_intensity, const Window& window, R_xlen_t number_types) {
  double integral(0);
  const auto volume(window.volume());
  for(R_xlen_t i(0); i < number_types; ++i) {
    integral += dominating_intensity[i] * volume;
  }
  return integral / static_cast<double>(number_types);
}

template<typename Configuration, typename Chain, typename Model, typename Lambda>
std::pair<bool, Configuration> update_LU_and_check_coalescence(Chain& chain, const Model& model, const Lambda& dominating_intensity, R_xlen_t number_types) {
  Configuration points_not_in_L(chain.get_last_configuration()), L;  // U starts with the end value of Z, L is an empty point process
  chain.iterate_forward_in_time([&points_not_in_L, &L, &model, &dominating_intensity, number_types](auto&& point, bool is_birth, auto mark) {
    if(is_birth) {
      // TODO: Might be able to reorganise, take log() and get speed-up.
      const auto normalized_papangelou_L(model.compute_papangelou(point, number_types, L) / dominating_intensity[get_type(point)]);
      if(normalized_papangelou_L < 0 || normalized_papangelou_L > 1) {
        Rcpp::stop("Did not correctly normalize the Papangelou intensity");
      }
      if(normalized_papangelou_L > mark) {
        const auto normalized_papangelou_U(model.compute_papangelou(point, number_types, L, points_not_in_L) / dominating_intensity[get_type(point)]);
        if(normalized_papangelou_U < 0 || normalized_papangelou_U > 1) {
          Rcpp::stop("Did not correctly normalize the Papangelou intensity");
        }
        if(normalized_papangelou_U > mark) {
          add_point(L, std::forward<decltype(point)>(point));
        } else {
          add_point(points_not_in_L, std::forward<decltype(point)>(point));
        }
      }
    } else {
      if(!remove_point(points_not_in_L, point)) {
        remove_point(L, point);
      }
    }
  });
  return std::pair<bool, Configuration>(empty(points_not_in_L), L);
}

} // namespace detail

template<typename Configuration, typename Model, typename Window>
inline auto simulate_coupling_from_the_past(const Model& model, const Window& window, R_xlen_t number_types) {
  const auto dominating_intensity(model.get_dominating_intensity(window, number_types));
  // TODO: Catch std::bad_alloc exceptions that happen when the upper bound to the Papangelou intensity is too large.
  Backwards_Markov_chain<Configuration> Z(rppp_single<Configuration>(window, dominating_intensity, number_types));

  const auto integral_of_dominating_intensity(detail::get_integral_of_dominating_intensity(dominating_intensity, window, number_types));
  const auto T0(Z.extend_until_T0(integral_of_dominating_intensity, window, number_types));
  while(true) {
    R_CheckUserInterrupt();
    const auto coalescence(detail::update_LU_and_check_coalescence<Configuration>(Z, model, dominating_intensity, number_types));
    if(coalescence.first) {
      return coalescence.second;
    }
    Z.extend_backwards(T0, integral_of_dominating_intensity, window, number_types);
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST
