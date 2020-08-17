#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_wrapper.hpp"
#include "configuration/make_R_configuration.hpp"

#include "simulation/rppp_single.hpp"

#include "utility/construct_if_missing.hpp"
#include "utility/get_list_or_first_element.hpp"
#include "utility/get_number_types.hpp"
#include "utility/make_default_types.hpp"
#include "utility/window.hpp"

inline SEXP rppp_helper(const ppjsdm::Window& window, const Rcpp::NumericVector& lambda, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, R_xlen_t point_types) {
  Rcpp::List samples(nsim);

  for(R_xlen_t i(0); i < nsim; ++i) {
    const auto sample(ppjsdm::rppp_single<ppjsdm::Configuration_wrapper>(window, lambda, point_types));
    samples[i] = ppjsdm::make_R_configuration(sample, types);
  }
  return ppjsdm::get_list_or_first_element(samples, nsim == 1 && drop);
}

// [[Rcpp::export]]
SEXP rppp_cpp(SEXP window, SEXP lambda, R_xlen_t nsim, SEXP types, bool drop, Rcpp::NumericVector mark_range) {
  const auto number_types(ppjsdm::get_number_types_and_check_conformance(1, lambda, types));
  lambda = ppjsdm::construct_if_missing<Rcpp::NumericVector>(lambda, 1., number_types);
  types = ppjsdm::make_types(types, number_types, lambda);
  const auto cpp_window(ppjsdm::Window(window, mark_range));
  return rppp_helper(cpp_window, lambda, nsim, types, drop, number_types);
}
