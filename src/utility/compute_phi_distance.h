#ifndef INCLUDE_PPJSDM_PHI_DISTANCE
#define INCLUDE_PPJSDM_PHI_DISTANCE

#include "point_manipulation.h"

#include <cmath> // std::sqrt

namespace ppjsdm {
namespace varphi {

class Identity {
public:
  static double apply(double square_distance, int, int) {
    return std::sqrt(square_distance);
  }
};

class Inverse_square {
public:
  static constexpr double apply(double square_distance, int, int) {
    return 1. / square_distance;
  }
};

class Strauss {
public:
  explicit Strauss(Rcpp::NumericMatrix radius) noexcept: square_radius_(radius) {
    for(auto& entry: square_radius_) {
      entry *= entry;
    }
  }
  double apply(double square_distance, int i, int j) const {
    if(square_distance <= square_radius_(i, j)) {
      return 1.;
    } else {
      return 0.;
    }
  }

private:
  Rcpp::NumericMatrix square_radius_;
};

} // namespace varphi

template<typename Varphi, typename Point>
double compute_phi_distance(const Point& point1, const Point& point2, const Varphi& varphi) {
  const auto delta_x(get_x(point1) - get_x(point2));
  const auto delta_y(get_y(point1) - get_y(point2));
  const auto square_distance(delta_x * delta_x + delta_y * delta_y);
  return varphi.apply(square_distance, get_type(point1), get_type(point2));
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_PHI_DISTANCE
