#ifndef INCLUDE_PPJSDM_CONFIGURATION_UTILITIES
#define INCLUDE_PPJSDM_CONFIGURATION_UTILITIES

#include <Rcpp.h>

#include <tuple> // std::pair

// inline void add_to_configuration(Rcpp::S4 configuration, Rcpp::NumericVector location, R_xlen_t type) {
//   auto x{Rcpp::NumericVector(configuration.slot("x"))};
//   auto y{Rcpp::NumericVector(configuration.slot("y"))};
//   auto types{Rcpp::IntegerVector(configuration.slot("types"))};
//
//   x.push_back(location[0]);
//   y.push_back(location[1]);
//   types.push_back(type + 1);
//
//   configuration.update(Rcpp::wrap(configuration));
//
//   // TODO: Change this aweful hack
//   configuration.slot("x") = x;
//   configuration.slot("y") = y;
//   configuration.slot("types") = types;
// }
//
// inline std::pair<Rcpp::NumericVector, R_xlen_t> remove_random_point(Rcpp::S4 configuration) {
//   auto x{Rcpp::NumericVector(configuration.slot("x"))};
//   auto y{Rcpp::NumericVector(configuration.slot("y"))};
//   auto types{Rcpp::IntegerVector(configuration.slot("types"))};
//
//   const R_xlen_t index{Rcpp::sample(x.size(), 1, false, R_NilValue, false)[0]};
//   const Rcpp::NumericVector saved_location{x[index], y[index]};
//   const R_xlen_t saved_type{types[index] - 1};
//   x.erase(index);
//   y.erase(index);
//   types.erase(index);
//
//   // TODO: Change this aweful hack
//   configuration.slot("x") = x;
//   configuration.slot("y") = y;
//   configuration.slot("types") = types;
//
//   return std::make_pair(saved_location, saved_type);
// }

[[nodiscard]] inline Rcpp::S4 make_configuration(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::IntegerVector types_vector, Rcpp::CharacterVector types) {
  types_vector.attr("class") = "factor";
  types_vector.attr("levels") = types;

  Rcpp::S4 configuration{"Configuration"};
  configuration.slot("x") = x;
  configuration.slot("y") = y;
  configuration.slot("types") = types_vector;
  return configuration;
}

class Configuration_wrapper {
public:
  explicit Configuration_wrapper(Rcpp::S4 configuration):
    x_{Rcpp::NumericVector(configuration.slot("x"))},
    y_{Rcpp::NumericVector(configuration.slot("y"))},
    types_{Rcpp::NumericVector(configuration.slot("types"))} {}

  [[nodiscard]] Rcpp::NumericVector x() const {
    return x_;
  }

  [[nodiscard]] Rcpp::NumericVector y() const {
    return y_;
  }

  [[nodiscard]] Rcpp::IntegerVector types() const {
    return types_;
  }

  [[nodiscard]] R_xlen_t get_number_points() const {
    return x_.size();
  }

private:
  Rcpp::NumericVector x_;
  Rcpp::NumericVector y_;
  Rcpp::IntegerVector types_;
};

#endif // INCLUDE_PPJSDM_CONFIGURATION_UTILITIES
