#ifndef INCLUDE_SHORT_RANGE_POTENTIALS
#define INCLUDE_SHORT_RANGE_POTENTIALS

#include <Rcpp.h>

#include "potential.hpp"
#include "../../point/point_manipulation.hpp"
#include "../../point/square_distance.hpp"
#include "../../utility/lightweight_matrix.hpp"

#include <algorithm> // std::max
#include <cmath> // std::sqrt, std::exp, std::log
#include <memory> // std::shared_ptr

namespace ppjsdm {
namespace potentials {

// Note: Public inheritance in order to inherit member variables if they exist in Varphi.
template<template<class> class Varphi, typename FloatType>
class Short_range_potential: public Varphi<FloatType>, public Potential<FloatType> {
private:
  // This matrix might store some elaborate function of the radius, and therefore might need
  // to store in FloatType for extra precision.
  Lightweight_square_matrix<FloatType> matrix_;
  using size_t = typename decltype(matrix_)::size_type;
  using V = Varphi<FloatType>;
public:
  explicit Short_range_potential(Rcpp::NumericMatrix radius): matrix_(radius.nrow()) {
    const size_t dim(radius.nrow());
    if(static_cast<size_t>(radius.ncol()) != dim) {
      Rcpp::stop("The matrix is not a square matrix, as was expected.");
    }
    for(size_t i(0); i < static_cast<size_t>(dim); ++i) {
      for(size_t j(0); j < static_cast<size_t>(dim); ++j) {
        matrix_(i, j) = V::set(radius(i, j));
      }
    }
  }

  FloatType apply(FloatType normalized_square_distance, int i, int j) const override {
    return V::apply(normalized_square_distance, matrix_(i, j));
  }

  FloatType get_square_lower_endpoint(int, int) const override {
    return static_cast<FloatType>(0.);
  }

  bool is_nonincreasing_after_lower_endpoint() const override {
    return V::is_nonincreasing_after_lower_endpoint;
  }
  bool is_two_valued() const override {
    return V::is_two_valued;
  }
};

template<typename FloatType>
struct Bump_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = true;
  static constexpr bool is_two_valued = false;

  static FloatType set(double radius) {
    return -static_cast<FloatType>(radius) * std::log(static_cast<FloatType>(2));
  }

  static FloatType apply(FloatType square_distance, FloatType constant) {
    return static_cast<FloatType>(1.) - std::exp(constant / std::sqrt(square_distance));
  }
};

template<typename FloatType>
struct Square_bump_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = true;
  static constexpr bool is_two_valued = false;

  static FloatType set(double radius) {
    return -static_cast<FloatType>(radius) * static_cast<FloatType>(radius) * std::log(static_cast<FloatType>(2));
  }

  static FloatType apply(FloatType square_distance, FloatType constant) {
    return static_cast<FloatType>(1.) - std::exp(constant / square_distance);
  }
};

template<typename FloatType>
struct Square_exponential_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = true;
  static constexpr bool is_two_valued = false;

  static FloatType set(double radius) {
    if(radius > 0) {
      return -std::log(static_cast<FloatType>(2)) / (static_cast<FloatType>(radius) * static_cast<FloatType>(radius));
    } else {
      return -std::numeric_limits<FloatType>::infinity();
    }
  }

  static FloatType apply(FloatType square_distance, FloatType constant) {
    return std::exp(constant * square_distance);
  }
};

template<typename FloatType>
struct Exponential_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = true;
  static constexpr bool is_two_valued = false;

  static FloatType set(double radius) {
    if(radius > 0) {
      return -std::log(static_cast<FloatType>(2)) / static_cast<FloatType>(radius);
    } else {
      return -std::numeric_limits<FloatType>::infinity();
    }
  }

  static FloatType apply(FloatType square_distance, FloatType constant) {
    return std::exp(constant * std::sqrt(square_distance));
  }
};

template<typename FloatType>
struct Linear_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = true;
  static constexpr bool is_two_valued = false;

  static FloatType set(double radius) {
    if(radius > 0) {
      return -static_cast<FloatType>(1.) / static_cast<FloatType>(radius);
    } else {
      return -std::numeric_limits<FloatType>::infinity();
    }
  }

  static FloatType apply(FloatType square_distance, FloatType constant) {
    return std::max<FloatType>(0., static_cast<FloatType>(1.) + constant * std::sqrt(square_distance));
  }
};

template<typename FloatType>
struct Strauss_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = true;
  static constexpr bool is_two_valued = true;

  static FloatType set(double radius) {
    return static_cast<FloatType>(radius) * static_cast<FloatType>(radius);
  }

  static FloatType apply(FloatType square_distance, FloatType constant) {
    if(square_distance <= constant) {
      return static_cast<FloatType>(1.);
    } else {
      return static_cast<FloatType>(0.);
    }
  }
};

template<typename FloatType>
using Exponential = Short_range_potential<Exponential_implementation, FloatType>;

template<typename FloatType>
using Square_exponential = Short_range_potential<Square_exponential_implementation, FloatType>;

template<typename FloatType>
using Bump = Short_range_potential<Bump_implementation, FloatType>;

template<typename FloatType>
using Square_bump = Short_range_potential<Square_bump_implementation, FloatType>;

template<typename FloatType>
using Strauss = Short_range_potential<Strauss_implementation, FloatType>;

template<typename FloatType>
using Linear = Short_range_potential<Linear_implementation, FloatType>;

} // namespace potentials

const constexpr char* const short_range_models[] = {
  "exponential",
  "square_exponential",
  "bump",
  "square_bump",
  "Geyer",
  "linear"
};

template<typename FloatType>
inline std::shared_ptr<const Potential<FloatType>> make_short_range_object(Rcpp::CharacterVector model,
                                                                           Rcpp::NumericMatrix radius) {
  const auto model_string(model[0]);
  if(model_string == short_range_models[0]) {
    return std::make_shared<potentials::Exponential<FloatType>>(radius);
  } else if(model_string == short_range_models[1]) {
    return std::make_shared<potentials::Square_exponential<FloatType>>(radius);
  } else if(model_string == short_range_models[2]) {
    return std::make_shared<potentials::Bump<FloatType>>(radius);
  } else if(model_string == short_range_models[3]) {
    return std::make_shared<potentials::Square_bump<FloatType>>(radius);
  } else if(model_string == short_range_models[4]) {
    return std::make_shared<potentials::Strauss<FloatType>>(radius);
  } else if(model_string == short_range_models[5]) {
    return std::make_shared<potentials::Linear<FloatType>>(radius);
  } else {
    Rcpp::stop("Incorrect model entered. A call to show_short_range_models() will show you the available choices.\n");
  }
}

} // namespace ppjsdm

#endif // INCLUDE_SHORT_RANGE_POTENTIALS
