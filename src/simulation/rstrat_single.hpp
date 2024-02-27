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
        const auto sample_x(xmin + (filling_x + unif_rand()) * delta_x[type]);
        const auto sample_y(ymin + (filling_y + unif_rand()) * delta_y[type]);
        if(window.is_in(sample_x, sample_y)) {
          add_point(configuration, Marked_point(sample_x, sample_y, type, window.draw_mark()));
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
  std::vector<decltype(window.volume())> delta(n.size());
  for(typename decltype(delta)::size_type i(0); i < delta.size(); ++i) {
    const auto adjusted_n(std::floor(std::sqrt(static_cast<decltype(window.volume())>(n[i]))));
    delta[i] = std::sqrt(volume) / adjusted_n;
  }
  return rstratpp_single<Configuration>(window, delta, delta);
}

} // namespace ppjsdm

#endif // INCLUDE_RSTRAT_SINGLE
