#ifndef INCLUDE_PPJSDM_PHI_DISTANCE
#define INCLUDE_PPJSDM_PHI_DISTANCE

#include "point_manipulation.h"

#include <cmath> // std::sqrt
#include <vector> // std::vector

namespace ppjsdm {
namespace varphi {

class Identity {
public:
  static double apply(double square_distance, R_xlen_t, R_xlen_t) {
    return std::sqrt(square_distance);
  }
};

class Inverse_square {
public:
  static constexpr double apply(double square_distance, R_xlen_t, R_xlen_t) {
    return 1. / square_distance;
  }
};

class Strauss {
private:
  auto access_square_radii(R_xlen_t i, R_xlen_t j) const {
    return square_radii_[i * dim_ + j];
  }
  void set_square_radii(R_xlen_t i, R_xlen_t j, double r) {
    square_radii_[i * dim_ + j] = r * r;
  }
public:
  explicit Strauss(Rcpp::NumericMatrix radius): dim_(radius.ncol()), square_radii_(dim_ * dim_) {
    for(R_xlen_t i(0); i < dim_; ++i) {
      for(R_xlen_t j(0); j < dim_; ++j) {
        set_square_radii(i, j, radius(i, j));
      }
    }
  }
  double apply(double square_distance, R_xlen_t i, R_xlen_t j) const {
    if(square_distance <= access_square_radii(i, j)) {
      return 1.;
    } else {
      return 0.;
    }
  }
private:
  R_xlen_t dim_;
  std::vector<double> square_radii_;
};

} // namespace varphi

template<typename Varphi, typename Point>
inline double compute_phi_distance(const Point& point1, const Point& point2, const Varphi& varphi) {
  const auto delta_x(get_x(point1) - get_x(point2));
  const auto delta_y(get_y(point1) - get_y(point2));
  const auto square_distance(delta_x * delta_x + delta_y * delta_y);
  return varphi.apply(square_distance, get_type(point1), get_type(point2));
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_PHI_DISTANCE
