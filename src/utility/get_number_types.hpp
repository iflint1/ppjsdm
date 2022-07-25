#ifndef INCLUDE_PPJSDM_GET_NUMBER_TYPES
#define INCLUDE_PPJSDM_GET_NUMBER_TYPES

#include <Rcpp.h>
#include <Rinternals.h>

namespace ppjsdm {
namespace detail {

// Deduce from object its "length", i.e., the number of types it implies
inline auto get_length_with_check(SEXP x) {
  if(Rf_isMatrix(x)) {
    const auto n(Rf_nrows(x));
    if(n != Rf_ncols(x)) {
      Rcpp::stop("Found a non-square matrix in the arguments.");
    }
    return n;
  } else if(Rf_isNewList(x)) {
    const auto x_list(Rcpp::as<Rcpp::List>(x));
    const auto n(Rcpp::as<Rcpp::NumericMatrix>(x_list[0]).nrow());
    if(n != Rcpp::as<Rcpp::NumericMatrix>(x_list[0]).ncol()) {
      Rcpp::stop("First argument is not a square matrix.");
    }
    const auto length_x(x_list.size());
    for(decltype(x_list.size()) i(1); i < length_x; ++i) {
      if(n != Rcpp::as<Rcpp::NumericMatrix>(x_list[i]).nrow() ||
         n != Rcpp::as<Rcpp::NumericMatrix>(x_list[i]).ncol()) {
        Rcpp::stop("Found a matrix with incompatible size in the arguments.");
      }
    }
    return n;
  } else {
    return Rf_length(x);
  }
}

// End of loop
inline R_xlen_t get_number_types_helper(R_xlen_t default_number_types, R_xlen_t number_types) {
  return number_types == 0 ? default_number_types : number_types;
}

template<typename... Args>
inline R_xlen_t get_number_types_helper(R_xlen_t default_number_types, R_xlen_t number_types, SEXP x, Args... args) {
  if(Rf_isNull(x)) {
    return get_number_types_helper(default_number_types, number_types, args...);
  } else {
    const auto length_x(get_length_with_check(x));
    if(number_types != 0) {
      if(length_x != number_types) {
        Rcpp::stop("Two of the given arguments have incompatible sizes.");
      } else {
        return get_number_types_helper(default_number_types, number_types, args...);
      }
    } else {
      return get_number_types_helper(default_number_types, length_x, args...);
    }
  }
}

} // namespace detail

template<typename... Args>
inline auto get_number_types_and_check_conformance(R_xlen_t default_number_types, Args... args) {
  return detail::get_number_types_helper(default_number_types, 0, args...);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_GET_NUMBER_TYPES
