#ifndef INCLUDE_PPJSDM_PHI_DISTANCE
#define INCLUDE_PPJSDM_PHI_DISTANCE

#include <cmath> // std::sqrt, std::exp, std::log
#include <vector> // std::vector

namespace ppjsdm {
namespace varphi {

// class Identity {
// public:
//   static double apply(double square_distance, int, int) {
//     return std::sqrt(square_distance);
//   }
//
//   template<typename Window>
//   static double get_maximum(const Window& window) {
//     return std::sqrt(window.square_diameter());
//   }
// };
//
// class Square {
// public:
//   static double apply(double square_distance, int, int) {
//     return square_distance;
//   }
//
//   template<typename Window>
//   static double get_maximum(const Window& window) {
//     return window.square_diameter();
//   }
// };

class Square_exponential {
private:
  auto access_lambda(int i, int j) const {
    return lambda_[i * dim_ + j];
  }

  void set_lambda(int i, int j, double r) {
    lambda_[i * dim_ + j] = -std::log(2) / (r * r);
  }
public:
  explicit Square_exponential(Rcpp::NumericMatrix radius): dim_(radius.ncol()), lambda_(dim_ * dim_) {
    for(R_xlen_t i(0); i < dim_; ++i) {
      for(R_xlen_t j(0); j < dim_; ++j) {
        set_lambda(i, j, radius(i, j));
      }
    }
  }

  double apply(double square_distance, int i, int j) const {
    return std::exp(access_lambda(i, j) * square_distance);
  }

  template<typename Window>
  static double get_maximum(const Window&) {
    return 1.0;
  }
private:
  R_xlen_t dim_;
  std::vector<double> lambda_;
};

class Exponential {
private:
  auto access_lambda(int i, int j) const {
    return lambda_[i * dim_ + j];
  }

  void set_lambda(int i, int j, double r) {
    lambda_[i * dim_ + j] = -std::log(2) / r;
  }
public:
  explicit Exponential(Rcpp::NumericMatrix radius): dim_(radius.ncol()), lambda_(dim_ * dim_) {
    for(R_xlen_t i(0); i < dim_; ++i) {
      for(R_xlen_t j(0); j < dim_; ++j) {
        set_lambda(i, j, radius(i, j));
      }
    }
  }

  double apply(double square_distance, int i, int j) const {
    return std::exp(access_lambda(i, j) * std::sqrt(square_distance));
  }

  template<typename Window>
  static double get_maximum(const Window&) {
    return 1.0;
  }
private:
  R_xlen_t dim_;
  std::vector<double> lambda_;
};

class Strauss {
private:
  auto access_square_radii(int i, int j) const {
    return square_radii_[i * dim_ + j];
  }

  void set_square_radii(int i, int j, double r) {
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
      return 1.;
    } else {
      return 0.;
    }
  }

  template<typename Window>
  static double get_maximum(const Window&) {
    return 1.;
  }
private:
  R_xlen_t dim_;
  std::vector<double> square_radii_;
};

} // namespace varphi
} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_PHI_DISTANCE
