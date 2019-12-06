#ifndef INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
#define INCLUDE_PPJSDM_RBINOMIALPP_SINGLE

#include "sample_from_window.h"
#include "window_utilities.h"

#include <Rcpp.h>

template<Window WindowType>
[[nodiscard]] inline Rcpp::S4 rbinomialpp_single(const Window_wrapper<WindowType>& window, Rcpp::IntegerVector n, Rcpp::CharacterVector types, R_xlen_t point_types, R_xlen_t total_number) {
  auto x{Rcpp::NumericVector(Rcpp::no_init(total_number))};
  auto y{Rcpp::NumericVector(Rcpp::no_init(total_number))};
  auto factors{Rcpp::IntegerVector(Rcpp::no_init(total_number))};
  factors.attr("class") = "factor";
  factors.attr("levels") = types;

  int fill{0};
  for(R_xlen_t j{0}; j < point_types; ++j) {
    const auto points_to_add{n[j]};
    sample_from_window(window, x, y, factors, points_to_add, fill, j + 1);
    fill += points_to_add;
  }

  Rcpp::S4 configuration{"Configuration"};
  configuration.slot("x") = x;
  configuration.slot("y") = y;
  configuration.slot("types") = factors;
  return configuration;
}


#endif // INCLUDE_PPJSDM_RBINOMIALPP_SINGLE
