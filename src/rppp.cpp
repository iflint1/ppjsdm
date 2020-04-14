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
#include "utility/window.hpp"

template<typename Lambda>
inline SEXP rppp_helper(const ppjsdm::Window& window, const Lambda& lambda, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, R_xlen_t point_types) {
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
//' @param mark_range Range of additional marks to give to the points.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
SEXP rppp(SEXP window = R_NilValue, SEXP lambda = R_NilValue, R_xlen_t nsim = 1, SEXP types = R_NilValue, bool drop = true, Rcpp::NumericVector mark_range = Rcpp::NumericVector::create(1., 1.)) {
  const auto number_types(ppjsdm::get_number_types_and_check_conformance(lambda, types));
  lambda = ppjsdm::construct_if_missing<Rcpp::NumericVector>(lambda, 1., number_types);
  types = ppjsdm::make_types(types, number_types, lambda);
  const auto cpp_window(ppjsdm::Window(window, mark_range));
  return ppjsdm::call_on_list_or_vector(lambda, [number_types, &cpp_window, nsim, &types, drop](const auto& l) {
    return rppp_helper(cpp_window, l, nsim, types, drop, number_types);
  });
}
