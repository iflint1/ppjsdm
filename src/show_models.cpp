#include <Rcpp.h>

#include "saturated_model/saturated_model.hpp"

//' Show the authorised short range models.
//'
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
void show_short_range_models() {
  Rcpp::Rcout << "The authorised models are: \n";
  for(const auto& m: ppjsdm::short_range_models) {
    Rcpp::Rcout << m << '\n';
  }
}

//' Show the authorised medium range models.
//'
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
void show_medium_range_models() {
  Rcpp::Rcout << "The authorised models are: \n";
  for(const auto& m: ppjsdm::medium_range_models) {
    Rcpp::Rcout << m << '\n';
  }
}
