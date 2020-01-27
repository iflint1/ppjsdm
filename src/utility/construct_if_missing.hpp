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
  SEXP operator()(R_xlen_t number_types, U val) {
    // The constructor that default-constructs to a given value does not work on Windows with RTools,
    // so do the construction by hand.
    Type new_value(number_types);
    for(R_xlen_t i(0); i < number_types; ++i) {
      new_value[i] = treat_as_real(val);
    }
    return Rcpp::wrap(new_value);
  }
};

template<>
struct make_type<Rcpp::NumericMatrix> {
  template<typename U>
  SEXP operator()(R_xlen_t number_types, U val) {
    Rcpp::NumericMatrix new_value(number_types, number_types);
    for(R_xlen_t i(0); i < number_types; ++i) {
      for(R_xlen_t j(0); j < number_types; ++j) {
        new_value(i, j) = treat_as_real(val);
      }
    }
    return Rcpp::wrap(new_value);
  }
};

} // namespace detail

template<typename Type, typename U>
inline SEXP construct_if_missing(R_xlen_t number_types, SEXP x, U def) {
  if(Rf_isNull(x)) {
    return detail::make_type<Type>{}(number_types, def);
  } else {
    // TODO: C++17 if constexpr would be better
    if(std::is_same<Type, Rcpp::NumericMatrix>::value) {
      if(!Rf_isMatrix(x)) {
        return detail::make_type<Type>{}(number_types, x);
      }
    }
    return x;
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CONSTRUCT_IF_MISSING
