#ifndef INCLUDE_RSTRAT_SINGLE
#define INCLUDE_RSTRAT_SINGLE

#include <Rinternals.h>

#include "../configuration/configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"
#include "../utility/window.hpp"

#include <cmath> // std::sqrt

namespace ppjsdm {

template<typename Configuration, typename Vector1, typename Vector2>
inline auto rstratpp_single(const Window& window,
                            const Vector1& delta_x,
                            const Vector2& delta_y) {
  Configuration configuration{};

  const auto number_types(delta_x.size());

  const auto xmin(window.xmin());
  const auto xmax(window.xmax());
  const auto ymin(window.ymin());
  const auto ymax(window.ymax());

  using filling_t = decltype(delta_x.size());

  filling_t maximum_size(0);
  for(filling_t type(0); type < number_types; ++type) {
    const auto nx = static_cast<filling_t>((xmax - xmin) / delta_x[type]);
    const auto ny = static_cast<filling_t>((ymax - ymin) / delta_y[type]);
    maximum_size += nx * ny;
  }
  reserve_if_possible(configuration, maximum_size);

  for(filling_t type(0); type < number_types; ++type) {
    const auto nx = static_cast<filling_t>((xmax - xmin) / delta_x[type]);
    const auto ny = static_cast<filling_t>((ymax - ymin) / delta_y[type]);
    for(filling_t filling_x(0); filling_x < nx; ++filling_x) {
      for(filling_t filling_y(0); filling_y < ny; ++filling_y) {
        constexpr const auto epsilon = static_cast<decltype(unif_rand())>(0.000001);
        // Assume the window is convex. If one of the corners of the sampling quadrat is in the window,
        // then it means by sampling in the quadrat, a.s., we have >0 probability of sampling a point in the window.
        // Continue drawing points until one such point is drawn.
        if(window.is_in(xmin + filling_x * delta_x[type] + epsilon, ymin + filling_y * delta_y[type] + epsilon) ||
           window.is_in(xmin + filling_x * delta_x[type] + epsilon, ymin + (filling_y + 1) * delta_y[type] - epsilon) ||
           window.is_in(xmin + (filling_x + 1) * delta_x[type] - epsilon, ymin + filling_y * delta_y[type] + epsilon) ||
           window.is_in(xmin + (filling_x + 1) * delta_x[type] - epsilon, ymin + (filling_y + 1) * delta_y[type] - epsilon)) {
          while(true) {
            const auto sample_x(xmin + (filling_x + unif_rand()) * delta_x[type]);
            const auto sample_y(ymin + (filling_y + unif_rand()) * delta_y[type]);
            if(window.is_in(sample_x, sample_y)) {
              add_point(configuration, Marked_point(sample_x, sample_y, type, window.draw_mark()));
              break;
            }
          }
        }
      }
    }
  }

  return configuration;
}

template<typename Configuration, typename Vector>
inline auto rstratpp_single(const Window& window,
                            const Vector& n) {
  const auto volume(window.volume());
  const auto outer_volume((window.xmax() - window.xmin()) * (window.ymax() - window.ymin()));
  std::vector<decltype(window.volume())> delta_x(n.size());
  std::vector<decltype(window.volume())> delta_y(n.size());
  for(typename decltype(delta_x)::size_type i(0); i < delta_x.size(); ++i) {
    const auto adjusted_n(std::round(std::sqrt(static_cast<decltype(window.volume())>(n[i]) * static_cast<decltype(window.volume())>(outer_volume) / static_cast<decltype(window.volume())>(volume))));
    delta_x[i] = (window.xmax() - window.xmin()) / adjusted_n;
    delta_y[i] = (window.ymax() - window.ymin()) / adjusted_n;
  }
  return rstratpp_single<Configuration>(window, delta_x, delta_y);
}

} // namespace ppjsdm

#endif // INCLUDE_RSTRAT_SINGLE
