#ifndef INCLUDE_PPJSDM_CONFIGURATION_WRAPPER
#define INCLUDE_PPJSDM_CONFIGURATION_WRAPPER

#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration_manipulation.h"
#include "../point/point_manipulation.h"

namespace ppjsdm {

class Configuration_wrapper {
private:
  class Marked_point_reference {
  public:
    Marked_point_reference(double& x, double& y, int& type):
    x_(x), y_(y), type_(type) {}

    Marked_point_reference& operator=(const Marked_point& point) {
      x_ = get_x(point);
      y_ = get_y(point);
      type_ = get_type(point) + 1;
      return *this;
    }

  private:
    double& x_;
    double& y_;
    int& type_;
  };

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

  Marked_point_reference operator[](R_xlen_t index) {
    return Marked_point_reference(x_[index], y_[index], types_[index]);
  }

  auto emplace_back(double x, double y, int type) {
    x_.push_back(x);
    y_.push_back(y);
    types_.push_back(type + 1);
  }

  auto push_back(const Marked_point& point) {
    emplace_back(get_x(point), get_y(point), get_type(point));
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
  static inline auto remove_point_by_index(Configuration_wrapper& configuration, R_xlen_t index) {
    auto point(configuration[index]);
    configuration.erase(index);
    return point;
  }

  static inline auto remove_random_point(Configuration_wrapper& configuration) {
    const auto index(Rcpp::sample(size(configuration), 1, false, R_NilValue, false)[0]);
    return remove_point_by_index(configuration, index);
  }
};

} // namespace traits
} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CONFIGURATION_WRAPPER
