#include <Rcpp.h>
#include <Rmath.h>

#include "call_on_list_or_vector.h"
#include "get_list_or_first_element.h"
#include "rbinomialpp_single.h"
#include "make_default_types.h"
#include "resolve_defaults.h"
#include "window_utilities.h"

namespace ppjsdm {

template<typename S, typename T>
inline SEXP rppp_helper(const S& window, const T& lambda, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, R_xlen_t point_types) {
  Rcpp::List samples(nsim);
  Rcpp::IntegerVector number_points(Rcpp::no_init(point_types));

  for(R_xlen_t i(0); i < nsim; ++i) {
    R_xlen_t total_number(0);
    for(R_xlen_t j(0); j < point_types; ++j) {
      const auto points_to_add(R::rpois(window.volume() * static_cast<double>(lambda[j])));
      number_points[j] = points_to_add;
      total_number += points_to_add;
    }

    samples[i] = rbinomialpp_single(window, number_points, types, point_types, total_number);
  }
  return get_list_or_first_element(samples, nsim, drop);
}

} // namespace ppjsdm

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
  types = ppjsdm::make_default_types(types, lambda, point_types);
  return ppjsdm::call_on_wrapped_window(window, [point_types, &lambda, nsim, &types, drop](const auto& w) {
    return ppjsdm::call_on_list_or_vector(lambda, [point_types, &w, nsim, &types, drop](const auto& l) {
      return ppjsdm::rppp_helper(w, l, nsim, types, drop, point_types);
    });
  });
}
