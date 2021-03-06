#ifndef INCLUDE_SATURATED_MODEL
#define INCLUDE_SATURATED_MODEL

#include <Rcpp.h>

#include "potentials/medium_range_potentials.hpp"
#include "potentials/potential.hpp"
#include "potentials/short_range_potentials.hpp"

#include <memory> // std::shared_ptr

namespace ppjsdm {

template<typename FloatType = double>
class Saturated_model {
public:
  using float_t = FloatType;

  Saturated_model(Rcpp::CharacterVector model,
                  Rcpp::NumericMatrix radius,
                  unsigned long long int saturation):
    object_(make_short_range_object<FloatType>(model, radius)),
    saturation_(saturation) {}

  Saturated_model(Rcpp::CharacterVector model,
                  Rcpp::NumericMatrix medium_range,
                  Rcpp::NumericMatrix long_range,
                  unsigned long long int saturation):
  object_(make_medium_range_object<FloatType>(model, medium_range, long_range)),
  saturation_(saturation) {}

  bool is_nonincreasing_after_lower_endpoint() const {
    return object_->is_nonincreasing_after_lower_endpoint();
  }

  bool is_two_valued() const {
    return object_->is_two_valued();
  }

  FloatType apply(FloatType normalized_square_distance, int i, int j) const {
    return object_->apply(normalized_square_distance, i, j);
  }

  FloatType get_square_lower_endpoint(int i, int j) const {
    return object_->get_square_lower_endpoint(i, j);
  }

  auto get_saturation() const {
    return saturation_;
  }

  FloatType get_maximum() const {
    return static_cast<FloatType>(saturation_);
  }

private:
  std::shared_ptr<const Potential<FloatType>> object_;
  unsigned long long int saturation_;
};

} // namespace ppjsdm

#endif // INCLUDE_SATURATED_MODEL
