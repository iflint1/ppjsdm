#ifndef INCLUDE_PPJSDM_POINT_MANIPULATION
#define INCLUDE_PPJSDM_POINT_MANIPULATION

#include <tuple> // std::tuple, std::get

namespace ppjsdm {

using Marked_point = std::tuple<double, double, int>;

namespace traits {

template<typename Point>
struct point_manipulation_defaults {
  static inline auto get_x(const Point& point) {
    return std::get<0>(point);
  }
  static inline auto get_y(const Point& point) {
    return std::get<1>(point);
  }
  static inline auto get_type(const Point& point) {
    return std::get<2>(point);
  }
};

template<typename Point>
struct point_manipulation: public point_manipulation_defaults<Point> {};

} // namespace traits

template<typename Point>
inline auto get_x(const Point& point) {
  return traits::point_manipulation<Point>::get_x(point);
}

template<typename Point>
inline auto get_y(const Point& point) {
  return traits::point_manipulation<Point>::get_y(point);
}

template<typename Point>
inline auto get_type(const Point& point) {
  return traits::point_manipulation<Point>::get_type(point);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_POINT_MANIPULATION
