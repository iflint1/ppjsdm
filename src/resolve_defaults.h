#ifndef INCLUDE_PPJSDM_RESOLVE_DEFAULTS
#define INCLUDE_PPJSDM_RESOLVE_DEFAULTS

#include <Rcpp.h>
#include <Rinternals.h>

namespace ppjsdm {
namespace detail {

inline R_xlen_t get_number_types_helper(R_xlen_t number_types) {
  return number_types == 0 ? 1 : number_types;
}

template<typename... Args>
inline R_xlen_t get_number_types_helper(R_xlen_t number_types, SEXP x, Args&... args) {
  if(Rf_isNull(x)) {
    return get_number_types_helper(number_types, args...);
  } else {
    R_xlen_t length_x;
    if(Rf_isMatrix(x)) {
      length_x = Rf_nrows(x);
      if(length_x != Rf_ncols(x)) {
        Rcpp::stop("Found a non-square matrix in the arguments.");
      }
    } else {
      length_x = Rf_length(x);
    }
    if(number_types != 0 && length_x != number_types) {
      Rcpp::stop("Two of the given arguments have incompatible sizes");
    }
    return get_number_types_helper(length_x, args...);
  }
}

} // namespace detail

template<typename... Args>
inline auto get_number_types_and_check_conformance(Args&... args) {
  return detail::get_number_types_helper(0, args...);
}

template<typename Type, typename U>
inline SEXP construct_if_missing(R_xlen_t number_types, SEXP x, U def) {
  if(Rf_isNull(x)) {
    // The constructor that default-constructs to a given value does not work on Windows with RTools,
    // so do the construction by hand.
    Type new_value(number_types);
    for(R_xlen_t i(0); i < number_types; ++i) {
      new_value[i] = def;
    }
    return Rcpp::wrap(new_value);
  } else {
    return x;
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RESOLVE_DEFAULTS
