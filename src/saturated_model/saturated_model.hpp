#ifndef INCLUDE_SATURATED_MODEL
#define INCLUDE_SATURATED_MODEL

#include <Rcpp.h>

#include "potentials/medium_range_potentials.hpp"
#include "potentials/short_range_potentials.hpp"

#include <memory> // std::shared_ptr

namespace ppjsdm {

class Saturated_model {
public:
  Saturated_model(Rcpp::CharacterVector model,
                  Rcpp::NumericMatrix radius,
                  unsigned long long int saturation):
  object_(make_short_range_object(model, radius)),
  saturation_(saturation) {}

  Saturated_model(Rcpp::CharacterVector model,
                  Rcpp::NumericMatrix medium_range,
                  Rcpp::NumericMatrix long_range,
                  unsigned long long int saturation):
  object_(make_medium_range_object(model, medium_range, long_range)),
  saturation_(saturation) {}

  bool is_nonincreasing_after_lower_endpoint() const {
    return object_->is_nonincreasing_after_lower_endpoint();
  }

  bool is_two_valued() const {
    return object_->is_two_valued();
  }

  double apply(double normalized_square_distance, int i, int j) const {
    return object_->apply(normalized_square_distance, i, j);
  }

  double get_square_lower_endpoint(int i, int j) const {
    return object_->get_square_lower_endpoint(i, j);
  }

  auto get_saturation() const {
    return saturation_;
  }

  double get_maximum() const {
    return static_cast<double>(saturation_);
  }

private:
  std::shared_ptr<const Potential> object_;
  unsigned long long int saturation_;
};

namespace detail {

enum class dispersionMethod {two_values, nonincreasing_after_lower_endpoint, generic};

template<template<dispersionMethod> class T, typename... Args>
inline auto dispatch_model(const Saturated_model& model,
                           Args&&... args) {
  if(model.is_two_valued()) {
    return T<dispersionMethod::two_values>{}(model, std::forward<Args>(args)...);
  } else if(model.is_nonincreasing_after_lower_endpoint()) {
    return T<dispersionMethod::nonincreasing_after_lower_endpoint>{}(model, std::forward<Args>(args)...);
  } else {
    return T<dispersionMethod::generic>{}(model, std::forward<Args>(args)...);
  }
}

} // namespace detail
} // namespace ppjsdm

#endif // INCLUDE_SATURATED_MODEL
