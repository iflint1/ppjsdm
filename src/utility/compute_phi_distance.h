#ifndef INCLUDE_PPJSDM_PHI_DISTANCE
#define INCLUDE_PPJSDM_PHI_DISTANCE

#include "point_manipulation.h"

#include <cmath> // std::sqrt

namespace ppjsdm {
namespace varphi {

class Identity {
public:
  static double apply(double square_distance) {
    return std::sqrt(square_distance);
  }
};

class Inverse_square {
public:
  static constexpr double apply(double square_distance) {
    return 1. / square_distance;
  }
};

class Strauss {
public:
  explicit Strauss(double radius) noexcept: square_radius_(radius * radius) {}
  double apply(double square_distance) const {
    if(square_distance <= square_radius_) {
      return 1.;
    } else {
      return 0.;
    }
  }

private:
  double square_radius_;
};

} // namespace varphi

template<typename V, typename Point>
double compute_phi_distance(const Point& point1, const Point& point2, const V& varphi) {
  const auto delta_x(get_x(point1) - get_x(point2));
  const auto delta_y(get_y(point1) - get_y(point2));
  const auto square_distance(delta_x * delta_x + delta_y * delta_y);
  return varphi.apply(square_distance);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_PHI_DISTANCE
