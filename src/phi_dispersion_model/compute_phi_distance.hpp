#ifndef INCLUDE_PPJSDM_PHI_DISTANCE
#define INCLUDE_PPJSDM_PHI_DISTANCE

#include "../point/point_manipulation.hpp"

#include <cmath> // std::sqrt, STD::EXP
#include <utility> // std::forward
#include <vector> // std::vector

namespace ppjsdm {
namespace varphi {

class Identity {
public:
  static double apply(double square_distance, int, int) {
    return std::sqrt(square_distance);
  }

  template<typename Window>
  double get_maximum(const Window& window) const {
    return std::sqrt(window.square_diameter());
  }
};

class Square {
public:
  static double apply(double square_distance, int, int) {
    return square_distance;
  }

  template<typename Window>
  double get_maximum(const Window& window) const {
    return window.square_diameter();
  }
};

class Exponential {
public:
  static double apply(double square_distance, int, int) {
    return std::exp(std::sqrt(square_distance)) - 1.0;
  }

  template<typename Window>
  double get_maximum(const Window& window) const {
    return std::exp(std::sqrt(window.square_diameter())) - 1.0;
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

  double apply(double square_distance, int i, int j) const {
    if(square_distance <= access_square_radii(i, j)) {
      return 0.;
    } else {
      return 1.;
    }
  }

  template<typename Window>
  double get_maximum(const Window&) const {
    return 1.;
  }
private:
  R_xlen_t dim_;
  std::vector<double> square_radii_;
};

template<typename V>
class Varphi: public V {
public:
  template<typename... Args>
  Varphi(Args&&... args): V(std::forward<Args>(args)...) {}

  template<typename Point>
  double compute_phi_distance(const Point& point1, const Point& point2) const {
    const auto delta_x(get_x(point1) - get_x(point2));
    const auto delta_y(get_y(point1) - get_y(point2));
    const auto square_distance(delta_x * delta_x + delta_y * delta_y);
    return V::apply(square_distance, get_type(point1), get_type(point2));
  }

  double varphi(double square_distance, int i, int j) const {
    return V::apply(square_distance, i, j);
  }

  template<typename Window>
  double get_maximum(const Window& window) const {
    return V::get_maximum(window);
  }
};

} // namespace varphi
} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_PHI_DISTANCE
