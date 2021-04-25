#ifndef INCLUDE_COMPUTE_DISPERSION
#define INCLUDE_COMPUTE_DISPERSION

#include <Rcpp.h>

#include "compute_dispersion_implementation.hpp"
#include "saturated_model.hpp"

#include "../utility/for_each_container.hpp"

#include <utility> // std::forward
#include <vector> // std::vector

namespace ppjsdm {
namespace detail {

template<typename AbstractDispersion, typename Point, typename... Configurations>
inline auto generic_dispersion_computation(const Saturated_model& varphi,
                                           const Point& point,
                                           R_xlen_t number_types,
                                           Configurations&... configurations) {
  // TODO: Something important to think about:
  // When running gibbsm or vcov.gibbsm, we need the saturation of every point in the configuration w.r.t.
  // each other point in the configuration.
  // In the current way I'm running this, we end up computing saturation of each point multiple times.
  // Large speed gains are possible by finding a better strategy.

  // TODO: The max size of the heap I'm using is always the saturation; reserve.
  using ValueType = typename AbstractDispersion::ValueType;
  using CountType = std::vector<ValueType>;
  using DispersionType = std::vector<double>;

  CountType count_vector(number_types);
  DispersionType dispersion(number_types);
  if(static_cast<decltype(size(configurations...))>(varphi.get_saturation()) >= size(configurations...)) {
    for_each_container([&count_vector, &point, &varphi](const auto& current_point) {
      if(!is_equal(current_point, point)) {
        AbstractDispersion::template update_count<std::numeric_limits<int>::max()>(varphi, count_vector[get_type(current_point)], current_point, point);
      }
    }, configurations...);
    add_count_to_dispersion<0, AbstractDispersion, 2>(varphi, dispersion, count_vector, point);
  } else {
    for_each_container([&dispersion, &count_vector, &point, &varphi,
                       saturation = varphi.get_saturation(), &configurations...](const auto& current_point) {
                         if(!is_equal(current_point, point)) {
                           ValueType count{};
                           // TODO: In the Geyer case, there's potential to use conditional_for_each_container and break if count == saturation.
                           // However, note that this would make the function different from other variations.
                           for_each_container([&point, &count, &current_point, &varphi, saturation](const auto& other_point) {
                             if(!is_equal(other_point, point) && !is_equal(other_point, current_point) && get_type(other_point) == get_type(point)) {
                               AbstractDispersion::template update_count<0>(varphi, count, current_point, other_point);
                             }
                           }, configurations...);
                           AbstractDispersion::template update_count<0>(varphi, count_vector[get_type(current_point)],
                                                                        current_point, point);
                           dispersion[get_type(current_point)] += AbstractDispersion::template delta<0, false>(varphi, count, current_point, point);
                         }
                       }, configurations...);
    add_count_to_dispersion<0, AbstractDispersion, 1>(varphi, dispersion, count_vector, point);
  }
  return dispersion;
}

template<detail::dispersionMethod T>
struct Functor_generic {
  template<typename... Args>
  auto operator()(Args&&... args) const {
    return generic_dispersion_computation<detail::compute_dispersion_implementation<T>>(std::forward<Args>(args)...);
  }
};

} // namespace detail

template<typename Point, typename... Configurations>
inline auto compute_dispersion(const Saturated_model& model,
                               const Point& point,
                               R_xlen_t number_types,
                               Configurations&&... configurations) {
  return detail::dispatch_model<detail::Functor_generic>(model, point, number_types, std::forward<Configurations>(configurations)...);
}

} // namespace ppjsdm

#endif // INCLUDE_COMPUTE_DISPERSION
