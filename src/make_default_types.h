#ifndef INCLUDE_PPJSDM_MAKE_CONFIGURATION
#define INCLUDE_PPJSDM_MAKE_CONFIGURATION

#include <Rcpp.h>

#include <string> // std::to_string

inline Rcpp::CharacterVector make_default_types(R_xlen_t size) {
  Rcpp::CharacterVector default_types(Rcpp::no_init(size));

  for(R_xlen_t i(0); i < size; ++i) {
    default_types[i] = std::to_string(i);
  }

  return default_types;
}


#endif // INCLUDE_PPJSDM_MAKE_CONFIGURATION
