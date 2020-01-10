#ifndef INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
#define INCLUDE_PPJSDM_RBINOMIALPP_SINGLE

#include <Rcpp.h>

#include "configuration_utilities.h"

namespace ppjsdm {

template<typename W, typename T>
inline auto rbinomialpp_single(const W& window, const T& n, Rcpp::CharacterVector types, R_xlen_t point_types, R_xlen_t total_number) {
  Rcpp::NumericVector x(Rcpp::no_init(total_number));
  Rcpp::NumericVector y(Rcpp::no_init(total_number));
  Rcpp::IntegerVector factors(Rcpp::no_init(total_number));

  R_xlen_t fill(0);
  for(R_xlen_t j(0); j < point_types; ++j) {
    const auto points_to_add(n[j]);
    for(R_xlen_t i(0); i < points_to_add; ++i) {
      const auto sample(window.sample());
      x[fill + i] = sample.first;
      y[fill + i] = sample.second;
      factors[fill + i] = j + 1;
    }
    fill += points_to_add;
  }

  return make_configuration(x, y, factors, types);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
