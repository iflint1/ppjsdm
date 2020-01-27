#ifndef INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST
#define INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST

#include <Rinternals.h>

#include <tuple> // std::pair

namespace ppjsdm {
namespace detail {

template<typename Lambda, typename Window>
inline auto get_integral_of_dominating_intensity(const Lambda& dominating_intensity, const Window& window) {
  Lambda copy(dominating_intensity);
  for(auto& i: copy) {
    i *= window.volume();
  }
  return copy;
}

template<typename Configuration, typename Chain, typename Model>
std::pair<bool, Configuration> update_LU_and_check_coalescence(Chain& chain, const Model& model, R_xlen_t point_types) {
  Configuration points_not_in_L(chain.get_last_configuration()), L;  // U starts with the end value of Z, L is an empty point process
  chain.iterate_forward_in_time([&points_not_in_L, &L](auto&& p, bool is_birth, auto mark) {
    if(is_birth) {
      const auto papangelou_L(model.compute_papangelou(L, p, point_types));
      if(papangelou > mark) {
        const auto papangelou_U(model.compute_papangelou(points_not_in_L, p, point_types) * papangelou_L);
        if(papangelou_U > mark) {
          add_point(L, std::forward<decltype(p)>(p));
        } else {
          add_point(points_not_in_L, std::forward<decltype(p)>(p));
        }
      }
    } else {
      if(!remove_point(points_not_in_L, p)) {
        remove_point(L, p);
      }
    }
  });
  return std::pair<bool, Configuration>(empty(points_not_in_L), L);
}

} // namespace detail

template<typename Configuration, typename Model, typename Window, typename Lambda>
inline auto simulate_coupling_from_the_past(const Model& model, const Window& window, const Lambda& dominating_intensity, R_xlen_t point_types) {
  Backwards_Markov_chain<Configuration> Z(rppp_single<Configuration>(window, dominating_intensity, point_types));
  const auto T0(Z.extend_until_T0(detail::get_integral_of_dominating_intensity(dominating_intensity, window), window, point_types));
  while(true) {
    const auto coalescence(detail::update_LU_and_check_coalescence<Configuration>(Z, model, point_types));
    if(coalescence.first) {
      return coalescence.second;
    }
    Z.extend_backwards(T0, detail::get_integral_of_dominating_intensity(dominating_intensity, window), window, point_types);
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST
