#ifndef INCLUDE_PPJSDM_PHI_DISTANCE
#define INCLUDE_PPJSDM_PHI_DISTANCE

#include <Rinternals.h>

#include <cmath> // std::sqrt, std::exp, std::log
#include <limits> // std::numeric_limits
#include <type_traits> // std::enable_if
#include <vector> // std::vector

namespace ppjsdm {
namespace varphi {

// Note: Public inheritance in order to inherit parameter `nonzero_value` if it exists in Varphi.
template<typename Varphi>
class Generic_varphi: public Varphi {
private:
  auto access_lambda(int i, int j) const {
    return lambda_[i * dim_ + j];
  }

  void set_lambda(int i, int j, double r) {
    lambda_[i * dim_ + j] = Varphi::set(r);
  }
public:
  explicit Generic_varphi(Rcpp::NumericMatrix radius): dim_(radius.ncol()), lambda_(dim_ * dim_) {
    for(R_xlen_t i(0); i < dim_; ++i) {
      for(R_xlen_t j(0); j < dim_; ++j) {
        set_lambda(i, j, radius(i, j));
      }
    }
  }

  double apply(double square_distance, int i, int j) const {
    return Varphi::apply(square_distance, access_lambda(i, j));
  }

  template<typename Window>
  static double get_maximum(const Window&) {
    return 1.0;
  }
private:
  R_xlen_t dim_;
  std::vector<double> lambda_;
};

struct Bump_implementation {
  static double set(double radius) {
    return -radius * std::log(2);
  }

  static double apply(double square_distance, double constant) {
    return 1.0 - std::exp(constant / std::sqrt(square_distance));
  }
};

struct Square_bump_implementation {
  static double set(double radius) {
    return -radius * radius * std::log(2);
  }

  static double apply(double square_distance, double constant) {
    return 1.0 - std::exp(constant / square_distance);
  }
};

struct Square_exponential_implementation {
  static double set(double radius) {
    return -std::log(2) / (radius * radius);
  }

  static double apply(double square_distance, double constant) {
    return std::exp(constant * square_distance);
  }
};

struct Exponential_implementation {
  static double set(double radius) {
    return -std::log(2) / radius;
  }

  static double apply(double square_distance, double constant) {
    return std::exp(constant * std::sqrt(square_distance));
  }
};

struct Strauss_implementation {
  static constexpr double nonzero_value = 1.0;

  static double set(double radius) {
    return radius * radius;
  }

  static double apply(double square_distance, double constant) {
    if(square_distance <= constant) {
      return 1.;
    } else {
      return 0.;
    }
  }
};

using Bump = Generic_varphi<Bump_implementation>;
using Square_bump = Generic_varphi<Square_bump_implementation>;
using Exponential = Generic_varphi<Exponential_implementation>;
using Square_exponential = Generic_varphi<Square_exponential_implementation>;
using Strauss = Generic_varphi<Strauss_implementation>;

} // namespace varphi

template<typename T, typename = void>
struct has_nonzero_value: std::false_type {};

template<typename T>
struct has_nonzero_value<T, decltype(static_cast<void>(T::nonzero_value), void())>: std::true_type {};

template<typename T>
constexpr bool has_nonzero_value_v = has_nonzero_value<T>::value;

template<typename Varphi>
inline std::enable_if_t<!has_nonzero_value_v<Varphi>, double> get_nonzero_value() {
  return 0.;
}

template<typename Varphi>
inline std::enable_if_t<has_nonzero_value_v<Varphi>, double> get_nonzero_value() {
  return Varphi::nonzero_value;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_PHI_DISTANCE
