#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_wrapper.hpp"
#include "configuration/make_R_configuration.hpp"

#include "simulation/rppp_single.hpp"

#include "utility/call_on_list_or_vector.hpp"

#include "utility/construct_if_missing.hpp"
#include "utility/get_list_or_first_element.hpp"
#include "utility/get_number_types.hpp"
#include "utility/make_default_types.hpp"
#include "utility/window_utilities.hpp"

template<typename Window, typename Lambda>
inline SEXP rppp_helper(const Window& window, const Lambda& lambda, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, R_xlen_t point_types) {
  Rcpp::List samples(nsim);

  for(R_xlen_t i(0); i < nsim; ++i) {
    const auto sample(ppjsdm::rppp_single<ppjsdm::Configuration_wrapper>(window, lambda, point_types));
    samples[i] = ppjsdm::make_R_configuration(sample, types);
  }
  return ppjsdm::get_list_or_first_element(samples, nsim, drop);
}

//' Sample a Poisson point processes
//'
//' @param window Simulation window. Default is a Rectangle window [0, 1]^2.
//' @param lambda A vector representing the intensities of the multipoint Poisson point processes.
//' Default is a vector of same size as types, filled with ones.
//' @param nsim Number of samples to generate. Default is 1.
//' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
//' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration.
//' Default is TRUE.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
SEXP rppp(SEXP window = R_NilValue, SEXP lambda = R_NilValue, R_xlen_t nsim = 1, SEXP types = R_NilValue, bool drop = true) {
  const auto point_types(ppjsdm::get_number_types_and_check_conformance(lambda, types));
  lambda = ppjsdm::construct_if_missing<Rcpp::NumericVector>(point_types, lambda, 1.);
  types = ppjsdm::make_types(types, point_types, lambda);
  return ppjsdm::call_on_wrapped_window(window, [point_types, &lambda, nsim, &types, drop](const auto& w) {
    return ppjsdm::call_on_list_or_vector(lambda, [point_types, &w, nsim, &types, drop](const auto& l) {
      return rppp_helper(w, l, nsim, types, drop, point_types);
    });
  });
}
