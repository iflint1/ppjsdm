#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_wrapper.hpp"
#include "configuration/make_R_configuration.hpp"

#include "simulation/rbinomialpp_single.hpp"

#include "utility/call_on_list_or_vector.hpp"
#include "utility/construct_if_missing.hpp"
#include "utility/get_list_or_first_element.hpp"
#include "utility/get_number_types.hpp"
#include "utility/make_default_types.hpp"
#include "utility/sum.hpp"
#include "utility/window_utilities.hpp"

template<typename N>
inline SEXP rbinomialpp_helper(const ppjsdm::Window& window, const N& n, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, R_xlen_t number_types) {
  const auto total_number(ppjsdm::sum<R_xlen_t>(n, number_types));

  Rcpp::List samples(nsim);
  for(R_xlen_t i(0); i < nsim; ++i) {
    const auto sample(ppjsdm::rbinomialpp_single<ppjsdm::Configuration_wrapper>(window, n, number_types, total_number));
    samples[i] = ppjsdm::make_R_configuration(sample, types);
  }
  return ppjsdm::get_list_or_first_element(samples, nsim, drop);
}

//' Sample a binomial point processes
//'
//' @param window Simulation window. Default is a Rectangle window [0, 1]^2.
//' @param n A vector representing the number of points of each types of the multipoint binomial point processes.
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
SEXP rbinomialpp(SEXP window = R_NilValue, SEXP n = R_NilValue, R_xlen_t nsim = 1, SEXP types = R_NilValue, bool drop = true, Rcpp::NumericVector mark_range = Rcpp::NumericVector::create(1., 1.)) {
  const auto number_types(ppjsdm::get_number_types_and_check_conformance(n, types));
  n = ppjsdm::construct_if_missing<Rcpp::IntegerVector>(n, 1, number_types);
  types = ppjsdm::make_types(types, number_types, n);
  const auto cpp_window(ppjsdm::get_window_from_R_object(window, mark_range));
  return ppjsdm::call_on_list_or_vector(n, [&cpp_window, nsim, &types, drop, number_types](const auto& m) {
    return rbinomialpp_helper(cpp_window, m, nsim, types, drop, number_types);
  });
}
