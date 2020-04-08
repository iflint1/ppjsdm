#ifndef INCLUDE_PPJSDM_POINT_MANIPULATION
#define INCLUDE_PPJSDM_POINT_MANIPULATION

#include <tuple> // std::tuple, std::get

namespace ppjsdm {

using Marked_point = std::tuple<double, double, int, double>;

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
  static inline auto get_mark(const Point& point) {
    return std::get<3>(point);
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

template<typename Point>
inline auto get_mark(const Point& point) {
  return traits::point_manipulation<Point>::get_mark(point);
}

template<typename Point, typename Other>
inline bool is_equal(const Point& point, const Other& other) {
  return get_x(point) == get_x(other)
          && get_y(point) == get_y(other)
          && get_type(point) == get_type(other)
          && get_mark(point) == get_mark(other);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_POINT_MANIPULATION
