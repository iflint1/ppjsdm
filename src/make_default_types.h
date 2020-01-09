#ifndef INCLUDE_PPJSDM_MAKE_CONFIGURATION
#define INCLUDE_PPJSDM_MAKE_CONFIGURATION

#include <Rcpp.h>

#include <string> // std::string, std::to_string

namespace ppjsdm {

inline Rcpp::CharacterVector make_default_types(SEXP might_contain_names, R_xlen_t size) {
  const SEXP potential_names(RCPP_GET_NAMES(might_contain_names));
  const auto given_names = Rf_isNull(potential_names)?
                            Rcpp::CharacterVector(size):
                            Rcpp::as<Rcpp::CharacterVector>(potential_names);
  Rcpp::CharacterVector default_types(Rcpp::no_init(size));
  for(R_xlen_t i(0); i < size; ++i) {
    if(given_names[i] != "") {
      default_types[i] = given_names[i];
    } else {
      default_types[i] = std::string("type").append(std::to_string(i));
    }
  }

  return default_types;
}

template<typename T>
inline Rcpp::CharacterVector make_default_types(Rcpp::Nullable<Rcpp::CharacterVector> types, const T& might_contain_names, R_xlen_t size) {
  if(types.isNull()) {
    return make_default_types(might_contain_names, size);
  } else {
    return types.as();
  }
}

template<typename T>
inline Rcpp::CharacterVector make_default_types(Rcpp::Nullable<Rcpp::CharacterVector> types, const T& might_contain_names) {
  return make_default_types(types, might_contain_names, might_contain_names.size());
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_MAKE_CONFIGURATION
