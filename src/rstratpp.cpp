#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_wrapper.hpp"
#include "configuration/make_R_configuration.hpp"

#include "simulation/rstrat_single.hpp"

#include "utility/construct_if_missing.hpp"
#include "utility/get_list_or_first_element.hpp"
#include "utility/get_number_types.hpp"
#include "utility/make_default_types.hpp"
#include "utility/sum.hpp"
#include "utility/window.hpp"

inline SEXP rstratpp_helper(const ppjsdm::Window& window, const Rcpp::NumericVector& delta_x, const Rcpp::NumericVector& delta_y, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop) {
  Rcpp::List samples(nsim);
  for(R_xlen_t i(0); i < nsim; ++i) {
    const auto sample(ppjsdm::rstratpp_single<ppjsdm::Configuration_wrapper>(window, delta_x, delta_y));
    samples[i] = ppjsdm::make_R_configuration(sample, types);
  }
  return ppjsdm::get_list_or_first_element(samples, nsim == 1 && drop);
}

// [[Rcpp::export]]
SEXP rstratpp_cpp(SEXP window, SEXP delta_x, SEXP delta_y, R_xlen_t nsim, SEXP types, bool drop, Rcpp::NumericVector mark_range) {
  const auto number_types(ppjsdm::get_number_types_and_check_conformance(1, delta_x, delta_y, types));
  delta_x = ppjsdm::construct_if_missing<Rcpp::NumericVector>(delta_x, 1, number_types);
  if(Rf_isNull(delta_y)) {
    delta_y = delta_x;
  }
  types = ppjsdm::make_types(types, number_types, delta_x, delta_y);
  const auto cpp_window(ppjsdm::Window(window, mark_range));
  return rstratpp_helper(cpp_window, delta_x, delta_y, nsim, types, drop);
}
