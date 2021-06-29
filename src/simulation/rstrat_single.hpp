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
                            const Vector1& nx,
                            const Vector2& ny) {
  Configuration configuration{};

  // TODO: Reserve sum(mx * ny)?

  const auto number_types(nx.size());

  const auto x1(window.xmin());
  const auto x2(window.xmax());
  const auto y1(window.ymin());
  const auto y2(window.ymax());

  auto iterator(configuration.begin());
  using filling_t = std::remove_cv_t<std::remove_reference_t<decltype(nx[0])>>;
  for(decltype(nx.size()) type(0); type < number_types; ++type) {
    const auto width((x2 - x1) / nx[type]);
    const auto height((y2 - y1) / ny[type]);

    for(filling_t filling_x(0); filling_x < nx[type]; ++filling_x) {
      for(filling_t filling_y(0); filling_y < ny[type]; ++filling_y) {
        const auto sample_x(x1 + (filling_x + unif_rand()) * width);
        const auto sample_y(y1 + (filling_y + unif_rand()) * height);
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
