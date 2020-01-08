#ifndef INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
#define INCLUDE_PPJSDM_RBINOMIALPP_SINGLE

#include <Rcpp.h>

#include "configuration_utilities.h"
#include "sample_from_window.h"

namespace ppjsdm {

template<typename W>
inline Rcpp::List rbinomialpp_single(const W& window, Rcpp::IntegerVector n, Rcpp::CharacterVector types, R_xlen_t point_types, R_xlen_t total_number) {
  Rcpp::NumericVector x(Rcpp::no_init(total_number));
  Rcpp::NumericVector y(Rcpp::no_init(total_number));
  Rcpp::IntegerVector factors(Rcpp::no_init(total_number));

  int fill(0);
  for(R_xlen_t j(0); j < point_types; ++j) {
    const auto points_to_add(n[j]);
    sample_from_window(window, x, y, factors, points_to_add, fill, j + 1);
    fill += points_to_add;
  }

  return make_configuration(x, y, factors, types);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
