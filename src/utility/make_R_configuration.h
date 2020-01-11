#ifndef INCLUDE_PPJSDM_MAKE_R_CONFIGURATION
#define INCLUDE_PPJSDM_MAKE_R_CONFIGURATION

#include <Rcpp.h>

#include <type_traits> //std::remove_cv_t, std::remove_reference_t

namespace ppjsdm {

inline Rcpp::List make_R_configuration(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::IntegerVector types_vector, Rcpp::CharacterVector types) {
  types_vector.attr("class") = "factor";
  types_vector.attr("levels") = types;

  Rcpp::List configuration = Rcpp::List::create(Rcpp::Named("x") = x, Rcpp::Named("y") = y, Rcpp::Named("types") = types_vector);
  configuration.attr("class") = "Configuration";
  return configuration;
}

template<typename Configuration>
inline Rcpp::List make_R_configuration(const Configuration& configuration, Rcpp::CharacterVector types) {
  const auto configuration_size(size(configuration));
  using size_type = std::remove_cv_t<std::remove_reference_t<decltype(configuration_size)>>;
  Rcpp::NumericVector x(Rcpp::no_init(configuration_size));
  Rcpp::NumericVector y(Rcpp::no_init(configuration_size));
  Rcpp::IntegerVector types_vector(Rcpp::no_init(configuration_size));
  for(size_type i(0); i < configuration_size; ++i) {
    x[i] = configuration.x(i);
    y[i] = configuration.y(i);
    // TODO: This does not work for Configuration_wrapper which should be treated separately anyhow.
    types_vector[i] = configuration.types(i) + 1;
  }
  return make_R_configuration(x, y, types_vector, types);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_MAKE_R_CONFIGURATION
