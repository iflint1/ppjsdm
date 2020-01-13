#ifndef INCLUDE_PPJSDM_GET_NUMBER_TYPES
#define INCLUDE_PPJSDM_GET_NUMBER_TYPES

#include <Rinternals.h>

namespace ppjsdm {
namespace detail {

inline auto get_number_types_with_check(SEXP x) {
  if(Rf_isMatrix(x)) {
    const auto n(Rf_nrows(x));
    if(n != Rf_ncols(x)) {
      Rcpp::stop("Found a non-square matrix in the arguments.");
    }
    return n;
  } else {
    return Rf_length(x);
  }
}

inline R_xlen_t get_number_types_helper(R_xlen_t number_types) {
  return number_types == 0 ? 1 : number_types;
}

template<typename... Args>
inline R_xlen_t get_number_types_helper(R_xlen_t number_types, SEXP x, Args&... args) {
  if(Rf_isNull(x)) {
    return get_number_types_helper(number_types, args...);
  } else {
    const auto length_x(get_number_types_with_check(x));
    if(number_types != 0 && length_x != number_types) {
      Rcpp::stop("Two of the given arguments have incompatible sizes");
    }
    return get_number_types_helper(length_x > number_types ? length_x: number_types, args...);
  }
}

} // namespace detail

template<typename... Args>
inline auto get_number_types_and_check_conformance(Args&... args) {
  return detail::get_number_types_helper(0, args...);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_GET_NUMBER_TYPES
