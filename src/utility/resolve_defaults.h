#ifndef INCLUDE_PPJSDM_RESOLVE_DEFAULTS
#define INCLUDE_PPJSDM_RESOLVE_DEFAULTS

#include <Rcpp.h>
#include <Rinternals.h>

#include <type_traits> // std::is_same

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

template<typename... Args>
inline auto get_number_types_and_check_conformance(Args&... args) {
  return detail::get_number_types_helper(0, args...);
}

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

#endif // INCLUDE_PPJSDM_RESOLVE_DEFAULTS
