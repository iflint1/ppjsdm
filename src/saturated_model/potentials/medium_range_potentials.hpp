#ifndef INCLUDE_PPJSDM_MEDIUM_RANGE_POTENTIALS
#define INCLUDE_PPJSDM_MEDIUM_RANGE_POTENTIALS

#include <Rcpp.h>

#include "traits.hpp"
#include "../../point/point_manipulation.hpp"
#include "../../point/square_distance.hpp"
#include "../../utility/lightweight_matrix.hpp"

#include <cmath> // std::sqrt, std::exp, std::log
#include <type_traits> // std::enable_if, std::false_type, std::true_type

namespace ppjsdm {
namespace potentials {

// Note: Public inheritance in order to inherit member variables if they exist in Varphi.
template<typename Varphi>
class Medium_range_potential: public Varphi {
private:
  Lightweight_square_matrix<double> medium_;
  Lightweight_square_matrix<double> long_;
  using size_t = typename decltype(medium_)::size_type;
protected:
  Medium_range_potential(Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range):
  medium_(medium_range.nrow()),
  long_(medium_range.nrow()){
    const size_t medium_dim(medium_range.nrow());
    if(static_cast<size_t>(medium_range.ncol()) != medium_dim
         || static_cast<size_t>(long_range.ncol()) != medium_dim
         || static_cast<size_t>(long_range.nrow()) != medium_dim) {
      Rcpp::stop("One of the matrices does not have the right dimensions.");
    }
    for(size_t i(0); i < static_cast<size_t>(medium_dim); ++i) {
      for(size_t j(0); j < static_cast<size_t>(medium_dim); ++j) {
        const auto med(medium_range(i, j));
        const auto lng(long_range(i, j));
        medium_(i, j) = Varphi::set_lower(med, lng);
        long_(i, j) = Varphi::set_upper(med, lng);
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

  template<typename V = Varphi, std::enable_if_t<has_square_lower_endpoint_v<V>>* = nullptr>
  double get_square_lower_endpoint(size_t i, size_t j) const {
    return get_square_lower_endpoint(medium_range(i, j));
  }
};

struct Half_square_exponential_implementation {
  static constexpr bool is_nonincreasing = false;
  static constexpr bool is_nonincreasing_after_lower_endpoint = true;

  static double set_lower(double lower, double) {
    return lower;
  }

  static double set_upper(double lower, double upper) {
    const auto delta(upper - lower);
    return -std::log(2) / (delta * delta);
  }

  static double get_square_lower_endpoint(double lower) {
    return lower * lower;
  }

  static double apply(double square_distance, double lower, double upper) {
    if(square_distance >= lower * lower) {
      const auto distance(std::sqrt(square_distance) - lower);
      return std::exp(upper * distance * distance);
    } else {
      return 0.;
    }
  }
};

struct Two_sided_Strauss_implementation {
  static constexpr double nonzero_value = 1.0;

  static double set_lower(double lower, double) {
    return lower * lower;
  }

  static double set_upper(double, double upper) {
    return upper * upper;
  }

  static double apply(double square_distance, double lower, double upper) {
    if(square_distance >= lower && square_distance <= upper) {
      return 1.;
    } else {
      return 0.;
    }
  }
};

struct Two_sided_linear_implementation {
  static constexpr bool is_nonincreasing = false;

  static double set_lower(double lower, double) {
    return lower;
  }

  static double set_upper(double, double upper) {
    return upper;
  }

  static double apply(double square_distance, double lower, double upper) {
    const auto distance(std::sqrt(square_distance));
    if(2 * distance <= (lower + upper)) {
      if(distance <= lower) {
        return 0.;
      } else {
        return 2. / (upper - lower) * (distance - lower);
      }
    } else {
      if(distance >= upper) {
        return 0.;
      } else {
        return 2. / (upper - lower) * (upper - distance);
      }
    }
  }
};

using Medium_range_square_exponential = Medium_range_potential<Half_square_exponential_implementation>;
using Medium_range_Geyer = Medium_range_potential<Two_sided_Strauss_implementation>;
using Medium_range_linear = Medium_range_potential<Two_sided_linear_implementation>;

} // namespace potentials
} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_MEDIUM_RANGE_POTENTIALS
