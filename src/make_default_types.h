#ifndef INCLUDE_PPJSDM_MAKE_CONFIGURATION
#define INCLUDE_PPJSDM_MAKE_CONFIGURATION

#include <Rcpp.h>
#include <Rinternals.h>

#include <string> // std::string, std::to_string

namespace ppjsdm {

inline SEXP make_default_types(SEXP might_contain_names, R_xlen_t size) {
  const SEXP potential_names(RCPP_GET_NAMES(might_contain_names));
  const auto given_names = Rf_isNull(potential_names)?
                            Rcpp::CharacterVector(size):
                            Rcpp::as<Rcpp::CharacterVector>(potential_names);
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

template<typename T>
inline SEXP make_default_types(SEXP types, const T& might_contain_names, R_xlen_t size) {
  if(Rf_isNull(types)) {
    return make_default_types(might_contain_names, size);
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

template<typename T>
inline SEXP make_default_types(SEXP types, const T& might_contain_names) {
  return make_default_types(types, might_contain_names, might_contain_names.size());
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_MAKE_CONFIGURATION
