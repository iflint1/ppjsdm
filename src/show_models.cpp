#include <Rcpp.h>

#include "phi_dispersion_model/compute_phi_dispersion.h"

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
