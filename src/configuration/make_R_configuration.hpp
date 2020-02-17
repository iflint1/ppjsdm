#ifndef INCLUDE_PPJSDM_MAKE_R_CONFIGURATION
#define INCLUDE_PPJSDM_MAKE_R_CONFIGURATION

#include <Rcpp.h>

#include "../configuration/configuration_manipulation.hpp"
#include "../configuration/configuration_wrapper.hpp"

#include <algorithm> // std::copy_n

namespace ppjsdm {

inline Rcpp::List make_R_configuration(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::IntegerVector types_vector, Rcpp::CharacterVector types) {
  types_vector.attr("class") = "factor";
  types_vector.attr("levels") = types;

  auto configuration = Rcpp::List::create(Rcpp::Named("x") = x, Rcpp::Named("y") = y, Rcpp::Named("types") = types_vector);
  configuration.attr("class") = "Configuration";
  return configuration;
}

inline auto make_R_configuration(const Configuration_wrapper& configuration, Rcpp::CharacterVector types) {
  return make_R_configuration(configuration.x(), configuration.y(), configuration.types(), types);
}

template<typename Configuration>
inline auto make_R_configuration(const Configuration& configuration, Rcpp::CharacterVector types) {
  const auto configuration_size(size(configuration));
  Configuration_wrapper r_configuration(configuration_size);
  std::copy_n(configuration.begin(), configuration_size, r_configuration.begin());
  return make_R_configuration(r_configuration, types);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_MAKE_R_CONFIGURATION
