#ifndef INCLUDE_PPJSDM_MAKE_DEFAULT_TYPES
#define INCLUDE_PPJSDM_MAKE_DEFAULT_TYPES

#include <Rcpp.h>
#include <Rinternals.h>

#include <string> // std::string, std::to_string
#include <utility> // std::forward

namespace ppjsdm {
namespace detail {

inline Rcpp::CharacterVector get_best_names(R_xlen_t size) {
  return Rcpp::CharacterVector(size);
}

template<typename... Args>
inline Rcpp::CharacterVector get_best_names(R_xlen_t size, SEXP might_contain_names, Args&&... other) {
  const SEXP potential_names(RCPP_GET_NAMES(might_contain_names));
  if(Rf_isNull(potential_names)) {
    return get_best_names(size, std::forward<Args>(other)...);
  } else {
    return Rcpp::as<Rcpp::CharacterVector>(potential_names);
  }
}

template<typename... Args>
inline SEXP make_default_types(R_xlen_t size, Args&&... might_contain_names) {
  Rcpp::CharacterVector given_names(get_best_names(size, std::forward<Args>(might_contain_names)...));
  Rcpp::CharacterVector default_types(Rcpp::no_init(size));
  for(R_xlen_t i(0); i < size; ++i) {
    if(given_names[i] != "") {
      default_types[i] = given_names[i];
    } else {
      default_types[i] = std::string("type").append(std::to_string(i + 1));
    }
  }

  return Rcpp::wrap(default_types);
}

} // namespace detail

template<typename... Args>
inline SEXP make_types(SEXP types, R_xlen_t size, Args&&... might_contain_names) {
  if(Rf_isNull(types)) {
    return detail::make_default_types(size, std::forward<Args>(might_contain_names)...);
  } else if (Rf_isNewList(types)) {
    Rcpp::List list(types);
    const auto length_types(list.size());
    Rcpp::CharacterVector new_types(length_types);
    for(R_xlen_t i(0); i < length_types; ++i) {
      new_types[i] =  Rf_asChar(list[i]);
    }
    return Rcpp::wrap(new_types);
  } else {
    return types;
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_MAKE_DEFAULT_TYPES
