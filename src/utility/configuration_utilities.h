#ifndef INCLUDE_PPJSDM_CONFIGURATION_UTILITIES
#define INCLUDE_PPJSDM_CONFIGURATION_UTILITIES

#include <Rcpp.h>

#include <tuple> // std::tuple
#include <vector> // std::vector

namespace ppjsdm {

inline Rcpp::List make_configuration(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::IntegerVector types_vector, Rcpp::CharacterVector types) {
  types_vector.attr("class") = "factor";
  types_vector.attr("levels") = types;

  Rcpp::List configuration = Rcpp::List::create(Rcpp::Named("x") = x, Rcpp::Named("y") = y, Rcpp::Named("types") = types_vector);
  configuration.attr("class") = "Configuration";
  return configuration;
}

class Configuration_wrapper {
public:
  explicit Configuration_wrapper(Rcpp::List configuration):
    x_(Rcpp::as<Rcpp::NumericVector>(configuration["x"])),
    y_(Rcpp::as<Rcpp::NumericVector>(configuration["y"])),
    types_(Rcpp::as<Rcpp::IntegerVector>(configuration["types"])) {}

  double x(R_xlen_t index) const {
    return x_[index];
  }

  double y(R_xlen_t index) const {
    return y_[index];
  }

  int types(R_xlen_t index) const {
    return types_[index];
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

  R_xlen_t get_number_points() const {
    return x_.size();
  }

private:
  Rcpp::NumericVector x_;
  Rcpp::NumericVector y_;
  Rcpp::IntegerVector types_;
};

namespace detail {

using Marked_point = std::tuple<double, double, int>;

} // namespace detail

class StdVector_configuration_wrapper {
  using vector_type = std::vector<detail::Marked_point>;
  using const_iterator = typename vector_type::const_iterator;
public:
  explicit StdVector_configuration_wrapper():
    points_{} {}

  double x(R_xlen_t index) const {
    return std::get<0>(points_[index]);
  }

  double y(R_xlen_t index) const {
    return std::get<1>(points_[index]);
  }

  int types(R_xlen_t index) const {
    return std::get<2>(points_[index]);
  }

  auto emplace_back(double x, double y, int types) {
    return points_.emplace_back(x, y, types);
  }

  auto erase(const_iterator i) {
    return points_.erase(i);
  }

  R_xlen_t get_number_points() const {
    return points_.size();
  }

private:
  vector_type points_;
};

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CONFIGURATION_UTILITIES
