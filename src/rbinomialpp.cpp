#include <Rcpp.h>

#include "call_on_list_or_vector.h"
#include "get_list_or_first_element.h"
#include "rbinomialpp_single.h"
#include "make_default_types.h"
#include "sum.h"
#include "window_utilities.h"

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
//' @param window The window.
//' @param n A vector representing the number of points of each types of the multipoint binomial point processes.
//' @param nsim Number of samples to generate.
//' @param types Types of the points.
//' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
SEXP rbinomialpp(SEXP window, SEXP n = Rcpp::IntegerVector::create(1), R_xlen_t nsim = 1, Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue, bool drop = true) {
  return ppjsdm::call_on_wrapped_window(window, [&n, nsim, &types, drop](const auto& w) {
    return ppjsdm::call_on_list_or_vector(n, [&w, nsim, &types, drop](const auto& m) {
      const auto point_types(m.size());
      return ppjsdm::rbinomialpp_helper(w, m, nsim, ppjsdm::make_default_types(types, m, point_types), drop, point_types);
    });
  });
}
