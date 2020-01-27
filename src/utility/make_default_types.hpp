#ifndef INCLUDE_PPJSDM_MAKE_DEFAULT_TYPES
#define INCLUDE_PPJSDM_MAKE_DEFAULT_TYPES

#include <Rcpp.h>
#include <Rinternals.h>

#include <string> // std::string, std::to_string

namespace ppjsdm {
namespace detail {

inline Rcpp::CharacterVector get_best_names(R_xlen_t size) {
  return Rcpp::CharacterVector(size);
}

template<typename... Args>
inline Rcpp::CharacterVector get_best_names(R_xlen_t size, SEXP might_contain_names, Args... others) {
  const SEXP potential_names(RCPP_GET_NAMES(might_contain_names));
  if(Rf_isNull(potential_names)) {
    return get_best_names(size, others...);
  } else {
    return Rcpp::as<Rcpp::CharacterVector>(potential_names);
  }
}

template<typename... Args>
inline SEXP make_default_types(R_xlen_t size, Args... might_contain_names) {
  Rcpp::CharacterVector given_names(get_best_names(size, might_contain_names...));
  Rcpp::CharacterVector default_types(Rcpp::no_init(size));
  for(R_xlen_t i(0); i < size; ++i) {
    const auto current_name(given_names[i]);
    if(current_name != "") {
      default_types[i] = current_name;
    } else {
      default_types[i] = std::string("type").append(std::to_string(i + 1));
    }
  }

  return Rcpp::wrap(default_types);
}

} // namespace detail

template<typename... Args>
inline SEXP make_types(SEXP types, R_xlen_t size, Args... might_contain_names) {
  if(Rf_isNull(types)) {
    return detail::make_default_types(size, might_contain_names...);
  } else if (Rf_isNewList(types)) {
    const auto list(Rcpp::as<Rcpp::List>(types));
    const auto length_types(list.size());
    Rcpp::CharacterVector new_types(length_types);
    for(R_xlen_t i(0); i < length_types; ++i) {
      new_types[i] =  Rf_asChar(list[i]);
    }
    return Rcpp::wrap(new_types);
  } else if (Rf_isVector(types)) {
    return types;
  } else {
    Rcpp::stop("The type of `types` is not convertible to a character vector.");
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_MAKE_DEFAULT_TYPES
