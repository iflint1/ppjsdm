#ifndef INCLUDE_PPJSDM_RESOLVE_DEFAULTS
#define INCLUDE_PPJSDM_RESOLVE_DEFAULTS

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
  } else { // TODO: This doesn't work with matrices.
    const auto length_x(LENGTH(x));
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
inline SEXP default_construct_if_missing(R_xlen_t number_types, SEXP x, U def) {
  if(Rf_isNull(x)) {
    return Rcpp::wrap(Type(number_types, def));
  } else {
    return x;
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RESOLVE_DEFAULTS
