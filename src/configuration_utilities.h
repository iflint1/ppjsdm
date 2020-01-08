#ifndef INCLUDE_PPJSDM_CONFIGURATION_UTILITIES
#define INCLUDE_PPJSDM_CONFIGURATION_UTILITIES

#include <Rcpp.h>

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

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CONFIGURATION_UTILITIES
