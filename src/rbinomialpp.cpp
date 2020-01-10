#include <Rcpp.h>

#include "utility/call_on_list_or_vector.h"
#include "utility/get_list_or_first_element.h"
#include "utility/rbinomialpp_single.h"
#include "utility/make_default_types.h"
#include "utility/resolve_defaults.h"
#include "utility/sum.h"
#include "utility/window_utilities.h"

namespace ppjsdm {

template<typename W, typename T>
inline SEXP rbinomialpp_helper(const W& window, const T& n, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, R_xlen_t point_types) {
  const auto total_number(ppjsdm::sum<int>(n, point_types));

  Rcpp::List samples(nsim);
  for(R_xlen_t i(0); i < nsim; ++i) {
    samples[i] = rbinomialpp_single(window, n, types, point_types, total_number);
  }
  return get_list_or_first_element(samples, nsim, drop);
}

} // namespace ppjsdm

//' Sample a binomial point processes
//'
//' @param window Simulation window. Default is a Rectangle window [0, 1]^2.
//' @param n A vector representing the number of points of each types of the multipoint binomial point processes.
//' Default is a vector of same size as types, filled with ones.
//' @param nsim Number of samples to generate. Default is 1.
//' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
//' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration.
//' Default is TRUE.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
SEXP rbinomialpp(SEXP window = R_NilValue, SEXP n = R_NilValue, R_xlen_t nsim = 1, SEXP types = R_NilValue, bool drop = true) {
  const auto point_types(ppjsdm::get_number_types_and_check_conformance(n, types));
  n = ppjsdm::construct_if_missing<Rcpp::IntegerVector>(point_types, n, 1);
  types = ppjsdm::make_types(types, point_types, n);
  return ppjsdm::call_on_wrapped_window(window, [&n, nsim, &types, drop, point_types](const auto& w) {
    return ppjsdm::call_on_list_or_vector(n, [&w, nsim, &types, drop, point_types](const auto& m) {
      return ppjsdm::rbinomialpp_helper(w, m, nsim, types, drop, point_types);
    });
  });
}
