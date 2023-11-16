#include "utility/make_default_types.hpp"

// [[Rcpp::export]]
SEXP make_types(SEXP types, R_xlen_t size, SEXP might_contain_name) {
  return(ppjsdm::detail::make_types(types, size, might_contain_name));
}

// [[Rcpp::export]]
SEXP make_types2(SEXP types, R_xlen_t size, SEXP might_contain_name, SEXP might_contain_name2) {
  return(ppjsdm::detail::make_types(types, size, might_contain_name, might_contain_name2));
}

// [[Rcpp::export]]
SEXP make_default_types(R_xlen_t size) {
  return(ppjsdm::detail::make_default_types(size));
}
