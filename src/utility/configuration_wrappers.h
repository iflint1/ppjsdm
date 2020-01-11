#ifndef INCLUDE_PPJSDM_CONFIGURATION_WRAPPERS
#define INCLUDE_PPJSDM_CONFIGURATION_WRAPPERS

#include <Rcpp.h>

#include "configuration_manipulation.h"
#include "point_manipulation.h"

namespace ppjsdm {

class Configuration_wrapper {
public:
  explicit Configuration_wrapper(Rcpp::List configuration):
    x_(Rcpp::as<Rcpp::NumericVector>(configuration["x"])),
    y_(Rcpp::as<Rcpp::NumericVector>(configuration["y"])),
    types_(Rcpp::as<Rcpp::IntegerVector>(configuration["types"])) {}
  explicit Configuration_wrapper(R_xlen_t size):
    x_(Rcpp::no_init(size)),
    y_(Rcpp::no_init(size)),
    types_(Rcpp::no_init(size)) {}

  auto operator[](R_xlen_t index) const {
    return Marked_point(x_[index], y_[index], types_[index] - 1);
  }

  void erase(R_xlen_t index) {
    x_.erase(x_.begin() + index);
    y_.erase(y_.begin() + index);
    types_.erase(types_.begin() + index);
  }

  Rcpp::NumericVector x() const {
    return x_;
  }

  Rcpp::NumericVector y() const {
    return y_;
  }

  Rcpp::IntegerVector types() const {
    return types_;
  }

  auto emplace_back(double x, double y, int types) {
    x_.push_back(x);
    y_.push_back(y);
    types_.push_back(types);
  }

  auto push_back(const Marked_point& point) {
    x_.push_back(std::get<0>(point));
    y_.push_back(std::get<1>(point));
    types_.push_back(std::get<2>(point));
  }

  auto size() const {
    return x_.size();
  }

private:
  Rcpp::NumericVector x_;
  Rcpp::NumericVector y_;
  Rcpp::IntegerVector types_;
};

namespace traits {

template<>
struct configuration_manipulation<Configuration_wrapper>: public configuration_manipulation_defaults<Configuration_wrapper> {
  static inline auto remove_random_point(Configuration_wrapper& configuration) {
    const auto index(Rcpp::sample(size(configuration), 1, false, R_NilValue, false)[0]);
    auto point(configuration[index]);
    configuration.erase(index);
    return point;
  }
};

} // namespace traits
} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CONFIGURATION_WRAPPERS
