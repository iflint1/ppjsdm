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

template<typename Configuration>
inline Rcpp::List make_configuration(const Configuration& configuration, Rcpp::CharacterVector types) {
  const auto configuration_size(configuration.get_number_points());
  Rcpp::NumericVector x(Rcpp::no_init(configuration_size));
  Rcpp::NumericVector y(Rcpp::no_init(configuration_size));
  Rcpp::IntegerVector types_vector(Rcpp::no_init(configuration_size));
  for(R_xlen_t i(0); i < configuration_size; ++i) {
    x[i] = configuration.x(i);
    y[i] = configuration.y(i);
    types_vector[i] = configuration.types(i);
  }
  return make_configuration(x, y, types_vector, types);
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

  auto emplace_back(double x, double y, int types) {
    x_.push_back(x);
    y_.push_back(y);
    types_.push_back(types);
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
public:
  StdVector_configuration_wrapper():
    points_{} {}
  // explicit StdVector_configuration_wrapper(const Configuration_wrapper& other) {
  //   const auto other_size(other.get_number_points());
  //   points_ = vector_type(other_size);
  //     for(R_xlen_t i(0); i < other_size; ++i) {
  //       points_[i] = std::make_tuple(other.x(i), other.y(i), other.types(i));
  //     }
  //   }

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

  auto erase(R_xlen_t index) {
    return points_.erase(points_.begin() + index);
  }

  R_xlen_t get_number_points() const {
    return points_.size();
  }

private:
  vector_type points_;
};

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CONFIGURATION_UTILITIES
