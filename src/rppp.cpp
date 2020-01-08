#include <Rcpp.h>
#include <Rmath.h>

#include "call_on_list_or_vector.h"
#include "rbinomialpp_single.h"
#include "make_default_types.h"
#include "window_utilities.h"

namespace ppjsdm {

template<typename S, typename T>
inline SEXP rppp_helper2(const S& window, T lambda, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, R_xlen_t point_types) {
  const auto volume(window.volume());

  Rcpp::List samples(nsim);
  Rcpp::IntegerVector number_points(Rcpp::no_init(point_types));

  for(R_xlen_t i(0); i < nsim; ++i) {
    R_xlen_t total_number(0);
    for(R_xlen_t j(0); j < point_types; ++j) {
      const auto points_to_add(R::rpois(volume * static_cast<double>(lambda[j])));
      number_points[j] = points_to_add;
      total_number += points_to_add;
    }

    samples[i] = rbinomialpp_single(window, number_points, types, point_types, total_number);
  }
  if(nsim == 1 && drop == true) {
    return Rcpp::wrap(samples[0]);
  } else {
    return Rcpp::wrap(samples);
  }
}

} // namespace ppjsdm

//' Sample a Poisson point processes
//'
//' @param window The window.
//' @param lambda A vector representing the intensities of the multipoint Poisson point processes.
//' @param nsim Number of samples to generate.
//' @param types Types of the points.
//' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
SEXP rppp(SEXP window, SEXP lambda = Rcpp::NumericVector::create(1), R_xlen_t nsim = 1, Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue, bool drop = true) {
  return ppjsdm::call_on_wrapped_window(window, [&lambda, nsim, &types, drop](const auto& w) {
    return ppjsdm::call_on_list_or_vector(lambda, [&w, nsim, &types, drop](const auto& l) {
      const auto point_types(l.size());
      const auto types_vector(ppjsdm::make_default_types(types, point_types));
      return ppjsdm::rppp_helper2(w, l, nsim, types_vector, drop, point_types);
    });
  });
}
