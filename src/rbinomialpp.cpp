#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_wrapper.hpp"
#include "configuration/make_R_configuration.hpp"

#include "simulation/rbinomialpp_single.hpp"

#include "utility/construct_if_missing.hpp"
#include "utility/get_list_or_first_element.hpp"
#include "utility/get_number_types.hpp"
#include "utility/make_default_types.hpp"
#include "utility/window.hpp"

inline SEXP rbinomialpp_helper(const ppjsdm::Window& window,
                               const Rcpp::NumericVector& n,
                               R_xlen_t nsim,
                               Rcpp::CharacterVector types,
                               bool drop) {
  Rcpp::List samples(nsim);
  for(R_xlen_t i(0); i < nsim; ++i) {
    const auto sample(ppjsdm::rbinomialpp_single<ppjsdm::Configuration_wrapper>(window, n));
    samples[i] = ppjsdm::make_R_configuration(sample, types);
  }
  return ppjsdm::get_list_or_first_element(samples, nsim == 1 && drop);
}

// [[Rcpp::export]]
SEXP rbinomialpp_cpp(SEXP window, SEXP n, R_xlen_t nsim, SEXP types, bool drop, Rcpp::NumericVector mark_range) {
  const auto number_types(ppjsdm::get_number_types_and_check_conformance(1, n, types));
  n = ppjsdm::construct_if_missing<Rcpp::IntegerVector>(n, 1, number_types);
  types = ppjsdm::detail::make_types(types, number_types, n);
  const auto cpp_window(ppjsdm::Window(window, mark_range));
  return rbinomialpp_helper(cpp_window, n, nsim, types, drop);
}
