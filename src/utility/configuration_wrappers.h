#ifndef INCLUDE_PPJSDM_CONFIGURATION_WRAPPERS
#define INCLUDE_PPJSDM_CONFIGURATION_WRAPPERS

#include <Rcpp.h>

#include <tuple> // std::tuple
#include <vector> // std::vector

namespace ppjsdm {

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

  auto push_back(const std::tuple<double, double, int>& point) {
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

  auto begin() const {
    return points_.begin();
  }

  void emplace_back(double x, double y, int types) {
    points_.emplace_back(x, y, types);
  }

  template<typename Point>
  auto push_back(const Point& point) {
    return points_.push_back(point);
  }

  template<typename Iterator>
  auto erase(Iterator iterator) {
    return points_.erase(iterator);
  }

  auto size() const {
    return points_.size();
  }

private:
  vector_type points_;
};

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CONFIGURATION_WRAPPERS
