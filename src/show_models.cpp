#include <Rcpp.h>

#include "saturated_model/potentials/short_range_potentials.hpp"
#include "saturated_model/potentials/medium_range_potentials.hpp"
#include "utility/array_size.hpp"

//' Show the authorised short range models.
//'
//' @export
//' @useDynLib ppjsdm, .registration = TRUE
//' @import Rcpp
// [[Rcpp::export]]
Rcpp::CharacterVector show_short_range_models() {
  Rcpp::CharacterVector ret(ppjsdm::array_size(ppjsdm::short_range_models));
  for(R_xlen_t i(0); i < ret.size(); ++i) {
    ret[i] = ppjsdm::short_range_models[i];
  }
  return ret;
}

//' Show the authorised medium range models.
//'
//' @export
//' @useDynLib ppjsdm, .registration = TRUE
//' @import Rcpp
// [[Rcpp::export]]
Rcpp::CharacterVector show_medium_range_models() {
  Rcpp::CharacterVector ret(ppjsdm::array_size(ppjsdm::medium_range_models));
  for(R_xlen_t i(0); i < ret.size(); ++i) {
    ret[i] = ppjsdm::medium_range_models[i];
  }
  return ret;
}
