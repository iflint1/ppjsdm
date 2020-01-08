#ifndef INCLUDE_PPJSDM_PHI_DISTANCE
#define INCLUDE_PPJSDM_PHI_DISTANCE

#include <cmath> // std::sqrt

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
  explicit Strauss(double radius) noexcept: square_radius_{radius * radius} {}
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

}

template<typename V>
double compute_phi_distance(double x1, double y1, double x2, double y2, const V& varphi) {
  const auto delta_x(x2 - x1);
  const auto delta_y(y2 - y1);
  const auto square_distance(delta_x * delta_x + delta_y * delta_y);
  return varphi.apply(square_distance);
}

#endif // INCLUDE_PPJSDM_PHI_DISTANCE
