#ifndef INCLUDE_PPJSDM_SHORT_RANGE_POTENTIALS
#define INCLUDE_PPJSDM_SHORT_RANGE_POTENTIALS

#include <Rcpp.h>

#include "../../point/point_manipulation.hpp"
#include "../../point/square_distance.hpp"
#include "../../utility/lightweight_matrix.hpp"

#include <algorithm> // std::max
#include <cmath> // std::sqrt, std::exp, std::log

namespace ppjsdm {
namespace potentials {

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

  double apply(double normalized_square_distance, size_t i, size_t j) const {
    return Varphi::apply(normalized_square_distance, matrix_(i, j));
  }

  template<typename Point, typename Other>
  double apply(const Point& point, const Other& other) const {
    return apply(normalized_square_distance(point, other), get_type(point), get_type(other));
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

  static double apply(double square_distance, double constant) {
    return std::exp(constant * square_distance);
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

struct Linear_implementation {
  static constexpr bool is_nonincreasing = true;

  static double set(double radius) {
    return -1. / radius;
  }

  static double apply(double square_distance, double constant) {
    return std::max<double>(0., 1. + constant * std::sqrt(square_distance));
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

using Bump = Short_range_potential<Bump_implementation>;
using Square_bump = Short_range_potential<Square_bump_implementation>;
using Exponential = Short_range_potential<Exponential_implementation>;
using Square_exponential = Short_range_potential<Square_exponential_implementation>;
using Strauss = Short_range_potential<Strauss_implementation>;
using Linear = Short_range_potential<Linear_implementation>;

} // namespace potentials
} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SHORT_RANGE_POTENTIALS
