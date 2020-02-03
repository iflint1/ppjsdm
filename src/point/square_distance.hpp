#ifndef INCLUDE_PPJSDM_SQUARE_DISTANCE
#define INCLUDE_PPJSDM_SQUARE_DISTANCE

#include "point_manipulation.hpp"

namespace ppjsdm {

template<typename Point, typename Other>
inline auto square_distance(const Point& point, const Other& other) {
  const auto delta_x(get_x(point) - get_x(other));
  const auto delta_y(get_y(point) - get_y(other));
  return delta_x * delta_x + delta_y * delta_y;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SQUARE_DISTANCE
