#ifndef INCLUDE_COMPUTE_DISPERSION
#define INCLUDE_COMPUTE_DISPERSION

#include <Rcpp.h>

#include "compute_dispersion_implementation.hpp"
#include "saturated_model.hpp"

#include "../configuration/get_number_points.hpp"
#include "../utility/for_each_container.hpp"

#include <utility> // std::forward
#include <vector> // std::vector

namespace ppjsdm {
namespace detail {

template<typename AbstractDispersion>
struct generic_dispersion_computation {
  template<typename Point, typename FloatType, typename... Configurations>
  auto operator()(const Saturated_model<FloatType>& varphi,
                  const Point& point,
                  R_xlen_t number_types,
                  Configurations&... configurations) const {
    // TODO: The max size of the heap I'm using is always the saturation; reserve.
    using ValueType = typename AbstractDispersion::ValueType;
    using CountType = std::vector<ValueType>;
    using DispersionType = std::vector<FloatType>;

    CountType count_vector(number_types);
    DispersionType dispersion(number_types);

    // If the saturation parameter is sufficiently large, we are guaranteed to never saturate.
    // In that case, a more efficient algorithm is available.
    const auto max_points_by_type(get_number_points_in_most_numerous_type(configurations...));
    if(static_cast<decltype(max_points_by_type)>(varphi.get_saturation()) >= max_points_by_type + 1) {
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
};

} // namespace detail

template<typename Point, typename FloatType, typename... Configurations>
inline auto compute_dispersion(const Saturated_model<FloatType>& model,
                               const Point& point,
                               R_xlen_t number_types,
                               Configurations&&... configurations) {
  return detail::dispatch_model<detail::generic_dispersion_computation>(model, point, number_types, std::forward<Configurations>(configurations)...);
}

} // namespace ppjsdm

#endif // INCLUDE_COMPUTE_DISPERSION
