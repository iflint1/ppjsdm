#ifndef INCLUDE_PPJSDM_PHI_DISTANCE
#define INCLUDE_PPJSDM_PHI_DISTANCE

#include <Rcpp.h>
#include <Rinternals.h>

#include "../point/point_manipulation.hpp"
#include "../point/square_distance.hpp"
#include "../utility/lightweight_matrix.hpp"

#include <cmath> // std::sqrt, std::exp, std::log
#include <type_traits> // std::enable_if, std::false_type, std::true_type
#include <vector> // std::vector

namespace ppjsdm {
namespace varphi {

// Note: Public inheritance in order to inherit member variables if they exist in Varphi.
template<typename Varphi>
class Short_range_potential: public Varphi {
private:
  Lightweight_square_matrix<double> matrix_;
  using size_t = typename decltype(matrix_)::size_type;
protected:
  explicit Short_range_potential(Rcpp::NumericMatrix radius): matrix_(radius.nrow()) {
    const size_t dim(radius.nrow());
    if(static_cast<size_t>(radius.ncol()) != dim) {
      Rcpp::stop("The matrix is not a square matrix, as was expected.");
    }
    for(size_t i(0); i < static_cast<size_t>(dim); ++i) {
      for(size_t j(0); j < static_cast<size_t>(dim); ++j) {
        matrix_(i, j) = Varphi::set(radius(i, j));
      }
    }
  }

  double apply(double square_distance, size_t i, size_t j) const {
    return Varphi::apply(square_distance, matrix_(i, j));
  }

  template<typename Point, typename Other>
  double apply(const Point& point, const Other& other) const {
    return apply(square_distance(point, other), get_type(point), get_type(other));
  }

  template<typename Window>
  constexpr static double get_maximum(const Window&) {
    return 1.0;
  }
};

template<typename Varphi>
class Medium_range_potential: public Varphi {
private:
  Lightweight_square_matrix<double> medium_;
  Lightweight_square_matrix<double> long_;
  using size_t = typename decltype(medium_)::size_type;
protected:
  Medium_range_potential(Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range):
  medium_(medium_range.nrow()),
  long_(long_range.nrow()){
    const size_t medium_dim(medium_range.nrow());
    if(static_cast<size_t>(medium_range.ncol()) != medium_dim
         || static_cast<size_t>(long_range.ncol()) != medium_dim
         || static_cast<size_t>(long_range.nrow()) != medium_dim) {
      Rcpp::stop("One of the matrices does not have the right dimensions.");
    }
    for(size_t i(0); i < static_cast<size_t>(medium_dim); ++i) {
      for(size_t j(0); j < static_cast<size_t>(medium_dim); ++j) {
        medium_(i, j) = Varphi::set_lower(medium_range(i, j));
        long_(i, j) = Varphi::set_upper(long_range(i, j));
      }
    }
  }

  double apply(double square_distance, size_t i, size_t j) const {
    return Varphi::apply(square_distance, medium_(i, j), long_(i, j));
  }

  template<typename Point, typename Other>
  double apply(const Point& point, const Other& other) const {
    return apply(square_distance(point, other), get_type(point), get_type(other));
  }

  template<typename Window>
  constexpr static double get_maximum(const Window&) {
    return 1.0;
  }
};

struct Bump_implementation {
  static constexpr bool is_nonincreasing = true;

  static double set(double radius) {
    return -radius * std::log(2);
  }

  static double apply(double square_distance, double constant) {
    return 1.0 - std::exp(constant / std::sqrt(square_distance));
  }
};

struct Square_bump_implementation {
  static constexpr bool is_nonincreasing = true;

  static double set(double radius) {
    return -radius * radius * std::log(2);
  }

  static double apply(double square_distance, double constant) {
    return 1.0 - std::exp(constant / square_distance);
  }
};

struct Square_exponential_implementation {
  static constexpr bool is_nonincreasing = true;

  static double set(double radius) {
    return -std::log(2) / (radius * radius);
  }

  static double set_lower(double radius) {
    return radius * radius;
  }

  static double set_upper(double radius) {
    return set(radius);
  }

  static double apply(double square_distance, double constant) {
    return std::exp(constant * square_distance);
  }

  static double apply(double square_distance, double lower, double upper) {
    if(square_distance >= lower) {
      const auto distance(std::sqrt(square_distance) - lower);
      return std::exp(upper * distance * distance);
    } else {
      return 0.;
    }
  }
};

struct Exponential_implementation {
  static constexpr bool is_nonincreasing = true;

  static double set(double radius) {
    return -std::log(2) / radius;
  }

  static double apply(double square_distance, double constant) {
    return std::exp(constant * std::sqrt(square_distance));
  }
};

struct Strauss_implementation {
  static constexpr double nonzero_value = 1.0;
  static constexpr bool is_nonincreasing = true;

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

using Bump = Short_range_potential<Bump_implementation>;
using Square_bump = Short_range_potential<Square_bump_implementation>;
using Exponential = Short_range_potential<Exponential_implementation>;
using Square_exponential = Short_range_potential<Square_exponential_implementation>;
using Strauss = Short_range_potential<Strauss_implementation>;

using Medium_range_square_exponential = Medium_range_potential<Square_exponential_implementation>;

} // namespace varphi

template<typename T, typename = void>
struct has_nonzero_value: std::false_type {};

template<typename T>
struct has_nonzero_value<T, decltype(static_cast<void>(T::nonzero_value), void())>: std::true_type {};

template<typename T>
constexpr bool has_nonzero_value_v = has_nonzero_value<T>::value;

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_PHI_DISTANCE
