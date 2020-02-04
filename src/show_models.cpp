#include <Rcpp.h>

#include "saturated_varphi_model/saturated_varphi_model.hpp"

//' Show the authorised models.
//'
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
void show_models() {
  Rcpp::Rcout << "The authorised models are: \n";
  for(const auto& m: ppjsdm::models) {
    Rcpp::Rcout << m << '\n';
  }
}
