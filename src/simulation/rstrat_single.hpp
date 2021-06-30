#ifndef INCLUDE_RSTRAT_SINGLE
#define INCLUDE_RSTRAT_SINGLE

#include <Rinternals.h>

#include <type_traits> // std::remove_cv_t, std::remove_reference_t

#include "../configuration/configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"
#include "../utility/window.hpp"

namespace ppjsdm {

template<typename Configuration, typename Vector1, typename Vector2>
inline auto rstratpp_single(const Window& window,
                            const Vector1& delta_x,
                            const Vector2& delta_y) {
  Configuration configuration{};

  // TODO: Reserve sum(mx * ny)?

  const auto number_types(delta_x.size());

  const auto x1(window.xmin());
  const auto x2(window.xmax());
  const auto y1(window.ymin());
  const auto y2(window.ymax());

  auto iterator(configuration.begin());
  using filling_t = decltype(delta_x.size());
  for(decltype(delta_x.size()) type(0); type < number_types; ++type) {
    const auto nx = static_cast<filling_t>((x2 - x1) / delta_x[type]);
    const auto ny = static_cast<filling_t>((y2 - y1) / delta_y[type]);

    for(filling_t filling_x(0); filling_x < nx; ++filling_x) {
      for(filling_t filling_y(0); filling_y < ny; ++filling_y) {
        const auto sample_x(x1 + (filling_x + unif_rand()) * delta_x[type]);
        const auto sample_y(y1 + (filling_y + unif_rand()) * delta_y[type]);
        if(window.is_in(sample_x, sample_y)) {
          add_point(configuration, Marked_point(sample_x, sample_y, type, window.draw_mark()));
        }
      }
    }
  }

  return configuration;
}

} // namespace ppjsdm

#endif // INCLUDE_RSTRAT_SINGLE
