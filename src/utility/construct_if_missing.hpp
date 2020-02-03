#ifndef INCLUDE_PPJSDM_CONSTRUCT_IF_MISSING
#define INCLUDE_PPJSDM_CONSTRUCT_IF_MISSING

#include <Rcpp.h>
#include <Rinternals.h>

#include <type_traits> // std::is_same

namespace ppjsdm {
namespace detail {

template<typename U>
inline auto treat_as_real(U value) {
  return value;
}

inline auto treat_as_real(SEXP value) {
  return Rf_asReal(value);
}

template<typename Type>
struct make_type {
  template<typename U>
  SEXP operator()(U val, R_xlen_t nrows, R_xlen_t) {
    // The constructor that default-constructs to a given value does not work on Windows with RTools,
    // so do the construction by hand.
    Type new_value(nrows);
    for(R_xlen_t i(0); i < nrows; ++i) {
      new_value[i] = treat_as_real(val);
    }
    return Rcpp::wrap(new_value);
  }
};

template<>
struct make_type<Rcpp::NumericMatrix> {
  template<typename U>
  SEXP operator()(U val, R_xlen_t nrows, R_xlen_t ncols) {
    Rcpp::NumericMatrix new_value(nrows, ncols);
    for(R_xlen_t i(0); i < nrows; ++i) {
      for(R_xlen_t j(0); j < ncols; ++j) {
        new_value(i, j) = treat_as_real(val);
      }
    }
    return Rcpp::wrap(new_value);
  }
};

} // namespace detail

template<typename Type, typename U>
inline SEXP construct_if_missing(SEXP x, U def, R_xlen_t nrows, R_xlen_t ncols) {
  if(Rf_isNull(x)) {
    return detail::make_type<Type>{}(def, nrows, ncols);
  } else {
    // TODO: C++17 if constexpr would be better
    if(std::is_same<Type, Rcpp::NumericMatrix>::value) {
      if(!Rf_isMatrix(x)) {
        return detail::make_type<Type>{}(x, nrows, ncols);
      }
    }
    return x;
  }
}

template<typename Type, typename U>
inline SEXP construct_if_missing(SEXP x, U def, R_xlen_t nrows) {
  return construct_if_missing<Type>(x, def, nrows, nrows);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CONSTRUCT_IF_MISSING
