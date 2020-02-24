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

template<typename Point, typename Other>
inline auto normalized_square_distance(const Point& point, const Other& other) {
  const auto delta_x(get_x(point) - get_x(other));
  const auto delta_y(get_y(point) - get_y(other));
  const auto sum_marks(get_mark(point) + get_mark(other));
  return (delta_x * delta_x + delta_y * delta_y) * 4.0 / (sum_marks * sum_marks);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SQUARE_DISTANCE
