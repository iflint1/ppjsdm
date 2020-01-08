#include <Rcpp.h>

#include "rbinomialpp_single.h"
#include "make_default_types.h"
#include "window_utilities.h"

template<typename W>
[[nodiscard]] inline SEXP rbinomialpp_helper(const W& window, Rcpp::IntegerVector n, R_xlen_t nsim, Rcpp::Nullable<Rcpp::CharacterVector> types, bool drop) {
  const auto point_types(n.size());
  const auto total_number(sum(n));

  if(types.isNull()) {
    types = Rcpp::wrap(make_default_types(point_types));
  }

  Rcpp::List samples(nsim);
  for(R_xlen_t i(0); i < nsim; ++i) {
    samples[i] = rbinomialpp_single(window, n, types.as(), point_types, total_number);
  }
  if(nsim == 1 && drop == true) {
    return Rcpp::wrap(samples[0]);
  } else {
    return Rcpp::wrap(samples);
  }
}


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
[[nodiscard]] SEXP rbinomialpp(SEXP window, Rcpp::IntegerVector n = Rcpp::IntegerVector::create(1), R_xlen_t nsim = 1, Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue, bool drop = true) {
  return call_on_wrapped_window(window, [&n, nsim, &types, drop](const auto& w) { return rbinomialpp_helper(w, n, nsim, types, drop); });
}
