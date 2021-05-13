#ifndef INCLUDE_MEDIUM_RANGE_POTENTIALS
#define INCLUDE_MEDIUM_RANGE_POTENTIALS

#include <Rcpp.h>

#include "potential.hpp"
#include "../../point/point_manipulation.hpp"
#include "../../point/square_distance.hpp"
#include "../../utility/lightweight_matrix.hpp"
#include "../../utility/window.hpp"

#include <cmath> // std::sqrt, std::exp, std::log
#include <memory> // std::shared_ptr

namespace ppjsdm {
namespace potentials {

// Note: Public inheritance in order to inherit member variables if they exist in Varphi.
template<template<class> class Varphi, typename FloatType>
class Medium_range_potential: public Varphi<FloatType>, public Potential<FloatType> {
private:
  // These matrices might store some elaborate function of the radii, and therefore might need
  // to store in FloatType for extra precision.
  Lightweight_square_matrix<FloatType> medium_;
  Lightweight_square_matrix<FloatType> long_;
  using size_t = typename decltype(medium_)::size_type;
  using V = Varphi<FloatType>;
public:
  Medium_range_potential(Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range):
  medium_(medium_range.nrow()),
  long_(medium_range.nrow()) {
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
        medium_(i, j) = V::set_lower(med, lng);
        long_(i, j) = V::set_upper(med, lng);
      }
    }
  }

  FloatType apply(FloatType normalized_square_distance, int i, int j) const override {
    return V::apply(normalized_square_distance, medium_(i, j), long_(i, j));
  }

  FloatType get_square_lower_endpoint(int i, int j) const override {
    return V::get_square_lower_endpoint(medium_(i, j));
  }

  bool is_nonincreasing_after_lower_endpoint() const override {
    return V::is_nonincreasing_after_lower_endpoint;
  }
  bool is_two_valued() const override {
    return V::is_two_valued;
  }
};

template<typename FloatType>
struct Half_square_exponential_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = true;
  static constexpr bool is_two_valued = false;

  static FloatType set_lower(double lower, double) {
    return static_cast<FloatType>(lower);
  }

  static FloatType set_upper(double lower, double upper) {
    const auto delta = static_cast<FloatType>(upper) - static_cast<FloatType>(lower);
    if(delta > 0.) {
      return -std::log(static_cast<FloatType>(2)) / (delta * delta);
    } else {
      return -std::numeric_limits<FloatType>::infinity();
    }
  }

  static FloatType get_square_lower_endpoint(FloatType lower) {
    return static_cast<FloatType>(lower) * static_cast<FloatType>(lower);
  }

  static FloatType apply(FloatType square_distance, FloatType lower, FloatType upper) {
    if(square_distance >= lower * lower) {
      const auto distance(std::sqrt(square_distance) - lower);
      return std::exp(upper * distance * distance);
    } else {
      return static_cast<FloatType>(0.);
    }
  }
};

template<typename FloatType>
struct Half_exponential_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = true;
  static constexpr bool is_two_valued = false;

  static FloatType set_lower(double lower, double) {
    return static_cast<FloatType>(lower);
  }

  static FloatType set_upper(double lower, double upper) {
    const auto delta = static_cast<FloatType>(upper) - static_cast<FloatType>(lower);
    if(delta > 0.) {
      return -std::log(static_cast<FloatType>(2)) / delta;
    } else {
      return -std::numeric_limits<FloatType>::infinity();
    }
  }

  static FloatType get_square_lower_endpoint(FloatType lower) {
    return static_cast<FloatType>(lower) * static_cast<FloatType>(lower);
  }

  static FloatType apply(FloatType square_distance, FloatType lower, FloatType upper) {
    if(square_distance >= lower * lower) {
      const auto distance(std::sqrt(square_distance) - lower);
      return std::exp(upper * distance);
    } else {
      return static_cast<FloatType>(0.);
    }
  }
};

template<typename FloatType>
struct Two_sided_square_exponential_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = false;
  static constexpr bool is_two_valued = false;

  static FloatType set_lower(double lower, double upper) {
    return (static_cast<FloatType>(lower) + static_cast<FloatType>(upper)) / static_cast<FloatType>(2.);
  }

  static FloatType set_upper(double lower, double upper) {
    const auto delta = static_cast<FloatType>(upper) - static_cast<FloatType>(lower);
    if(delta > 0) {
      return -static_cast<FloatType>(4) * std::log(static_cast<FloatType>(2)) / (delta * delta);
    } else {
      return -std::numeric_limits<FloatType>::infinity();
    }
  }

  static FloatType get_square_lower_endpoint(FloatType) {
    return static_cast<FloatType>(0.);
  }

  static FloatType apply(FloatType square_distance, FloatType lower, FloatType upper) {
    const auto distance(std::sqrt(square_distance) - lower);
    return std::exp(upper * distance * distance);
  }
};

template<typename FloatType>
struct Two_sided_exponential_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = false;
  static constexpr bool is_two_valued = false;

  static FloatType set_lower(double lower, double upper) {
    return (static_cast<FloatType>(lower) + static_cast<FloatType>(upper)) / static_cast<FloatType>(2.);
  }

  static FloatType set_upper(double lower, double upper) {
    const auto delta = static_cast<FloatType>(upper - lower);
    if(delta > 0.) {
      return -static_cast<FloatType>(2) * std::log(static_cast<FloatType>(2)) / delta;
    } else {
      return -std::numeric_limits<FloatType>::infinity();
    }
  }

  static FloatType get_square_lower_endpoint(FloatType) {
    return static_cast<FloatType>(0.);
  }

  static FloatType apply(FloatType square_distance, FloatType lower, FloatType upper) {
    const auto distance(std::abs(std::sqrt(square_distance) - lower));
    return std::exp(upper * distance);
  }
};

template<typename FloatType>
struct Two_sided_bump_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = false;
  static constexpr bool is_two_valued = false;

  static FloatType set_lower(double lower, double upper) {
    return (static_cast<FloatType>(lower) + static_cast<FloatType>(upper)) / static_cast<FloatType>(2.);
  }

  static FloatType set_upper(double lower, double upper) {
    return -static_cast<FloatType>(0.5) * std::log(static_cast<FloatType>(2)) * (static_cast<FloatType>(upper) - static_cast<FloatType>(lower));
  }

  static FloatType get_square_lower_endpoint(FloatType) {
    return static_cast<FloatType>(0.);
  }

  static FloatType apply(FloatType square_distance, FloatType lower, FloatType upper) {
    const auto distance(std::sqrt(square_distance));
    if(distance == lower) {
      return static_cast<FloatType>(1.);
    } else if(distance < lower) {
      return static_cast<FloatType>(1.) - std::exp(upper / (lower - distance));
    } else {
      return static_cast<FloatType>(1.) - std::exp(upper / (distance - lower));
    }
  }
};

template<typename FloatType>
struct Two_sided_square_bump_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = false;
  static constexpr bool is_two_valued = false;

  static FloatType set_lower(double lower, double upper) {
    return (static_cast<FloatType>(lower) + static_cast<FloatType>(upper)) / static_cast<FloatType>(2.);
  }

  static FloatType set_upper(double lower, double upper) {
    const auto delta = static_cast<FloatType>(upper) - static_cast<FloatType>(lower);
    return -static_cast<FloatType>(0.25) * std::log(static_cast<FloatType>(2)) * delta * delta;
  }

  static FloatType get_square_lower_endpoint(FloatType) {
    return static_cast<FloatType>(0.);
  }

  static FloatType apply(FloatType square_distance, FloatType lower, FloatType upper) {
    const auto distance(std::sqrt(square_distance));
    if(distance == lower) {
      return static_cast<FloatType>(1.);
    } else {
      const auto delta = static_cast<FloatType>(upper) - static_cast<FloatType>(lower);
      return static_cast<FloatType>(1.) - std::exp(upper / (delta * delta));
    }
  }
};

template<typename FloatType>
struct Two_sided_Strauss_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = false;
  static constexpr bool is_two_valued = true;

  static FloatType set_lower(double lower, double) {
    return static_cast<FloatType>(lower) * static_cast<FloatType>(lower);
  }

  static FloatType set_upper(double, double upper) {
    return static_cast<FloatType>(upper) * static_cast<FloatType>(upper);
  }

  static FloatType get_square_lower_endpoint(FloatType) {
    return static_cast<FloatType>(0.);
  }

  static FloatType apply(FloatType square_distance, FloatType lower, FloatType upper) {
    if(square_distance >= lower && square_distance <= upper) {
      return static_cast<FloatType>(1.);
    } else {
      return static_cast<FloatType>(0.);
    }
  }
};

template<typename FloatType>
struct Two_sided_linear_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = false;
  static constexpr bool is_two_valued = false;

  static FloatType set_lower(double lower, double) {
    return static_cast<FloatType>(lower);
  }

  static FloatType set_upper(double, double upper) {
    return static_cast<FloatType>(upper);
  }

  static FloatType get_square_lower_endpoint(FloatType) {
    return static_cast<FloatType>(0.);
  }

  static FloatType apply(FloatType square_distance, FloatType lower, FloatType upper) {
    const auto distance(std::sqrt(square_distance));
    if(2 * distance <= (lower + upper)) {
      if(distance <= lower || upper == lower) {
        return static_cast<FloatType>(0.);
      } else {
        return static_cast<FloatType>(2.) / (upper - lower) * (distance - lower);
      }
    } else {
      if(distance >= upper || upper == lower) {
        return static_cast<FloatType>(0.);
      } else {
        return static_cast<FloatType>(2.) / (upper - lower) * (upper - distance);
      }
    }
  }
};

template<typename FloatType>
struct Tanh_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = false;
  static constexpr bool is_two_valued = false;
  static constexpr const FloatType lambda = 5.;
  static constexpr const FloatType constant = 0.50678365490630417067308144397L; // Computed as 1 / (2 * tanh(lambda / 2))

  static FloatType set_lower(double lower, double) {
    return static_cast<FloatType>(lower);
  }

  static FloatType set_upper(double, double upper) {
    return static_cast<FloatType>(upper);
  }

  static FloatType get_square_lower_endpoint(FloatType) {
    return static_cast<FloatType>(0.);
  }

  static FloatType apply(FloatType square_distance, FloatType lower, FloatType upper) {
    const auto distance(std::sqrt(square_distance));
    const auto factor(lambda / (upper - lower));
    return constant * (std::tanh(factor * (distance - lower)) + std::tanh(factor * (upper - distance)));
  }
};

template<typename FloatType>
using Medium_range_square_exponential = Medium_range_potential<Two_sided_square_exponential_implementation, FloatType>;

template<typename FloatType>
using Medium_range_half_square_exponential = Medium_range_potential<Half_square_exponential_implementation, FloatType>;

template<typename FloatType>
using Medium_range_Geyer = Medium_range_potential<Two_sided_Strauss_implementation, FloatType>;

template<typename FloatType>
using Medium_range_linear = Medium_range_potential<Two_sided_linear_implementation, FloatType>;

template<typename FloatType>
using Medium_range_half_exponential = Medium_range_potential<Half_exponential_implementation, FloatType>;

template<typename FloatType>
using Medium_range_exponential = Medium_range_potential<Two_sided_exponential_implementation, FloatType>;

template<typename FloatType>
using Medium_range_bump = Medium_range_potential<Two_sided_bump_implementation, FloatType>;

template<typename FloatType>
using Medium_range_square_bump = Medium_range_potential<Two_sided_square_bump_implementation, FloatType>;

template<typename FloatType>
using Medium_range_tanh = Medium_range_potential<Tanh_implementation, FloatType>;

} // namespace potentials

const constexpr char* const medium_range_models[] = {
  "square_exponential",
  "half_square_exponential",
  "Geyer",
  "linear",
  "half_exponential",
  "exponential",
  "bump",
  "square_bump",
  "tanh"
};

template<typename FloatType>
inline std::shared_ptr<const Potential<FloatType>> make_medium_range_object(Rcpp::CharacterVector model,
                                                                            Rcpp::NumericMatrix medium_range,
                                                                            Rcpp::NumericMatrix long_range) {
  const auto model_string(model[0]);
  if(model_string == medium_range_models[0]) {
    return std::make_shared<potentials::Medium_range_square_exponential<FloatType>>(medium_range, long_range);
  } else if(model_string == medium_range_models[1]) {
    return std::make_shared<potentials::Medium_range_half_square_exponential<FloatType>>(medium_range, long_range);
  } else if(model_string == medium_range_models[2]) {
    return std::make_shared<potentials::Medium_range_Geyer<FloatType>>(medium_range, long_range);
  } else if(model_string == medium_range_models[3]) {
    return std::make_shared<potentials::Medium_range_linear<FloatType>>(medium_range, long_range);
  } else if(model_string == medium_range_models[4]) {
    return std::make_shared<potentials::Medium_range_half_exponential<FloatType>>(medium_range, long_range);
  } else if(model_string == medium_range_models[5]) {
    return std::make_shared<potentials::Medium_range_exponential<FloatType>>(medium_range, long_range);
  } else if(model_string == medium_range_models[6]) {
    return std::make_shared<potentials::Medium_range_bump<FloatType>>(medium_range, long_range);
  } else if(model_string == medium_range_models[7]) {
    return std::make_shared<potentials::Medium_range_square_bump<FloatType>>(medium_range, long_range);
  } else if(model_string == medium_range_models[8]) {
    return std::make_shared<potentials::Medium_range_tanh<FloatType>>(medium_range, long_range);
  } else {
    Rcpp::stop("Incorrect model entered. A call to show_medium_range_models() will show you the available choices.\n");
  }
}

} // namespace ppjsdm

#endif // INCLUDE_MEDIUM_RANGE_POTENTIALS
