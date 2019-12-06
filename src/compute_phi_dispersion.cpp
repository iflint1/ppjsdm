#include <Rcpp.h>

#include "compute_phi_dispersion.h"

//' Compute the phi-dispersion of a marked configuration (radius = 0 corresponds to the proposal model).
//'
//' @param configuration The configuration.
//' @param number_types Number of different types.
//' @param radius Radius of interaction.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
[[nodiscard]] Rcpp::NumericMatrix compute_phi_dispersion(Rcpp::S4 configuration, int number_types, double radius = 0) {
  if(radius == 0) {
    return compute_phi_dispersion_helper<Varphi::identity>(configuration, number_types);
  } else {
    return compute_phi_dispersion_helper<Varphi::Strauss>(configuration, number_types, radius * radius);
  }
}

//' Compute the delta of the phi-dispersion of a marked configuration (radius = 0 corresponds to the proposal model).
//'
//' @param configuration The configuration.
//' @param location Point to be added.
//' @param type Type of point to be added.
//' @param number_types Number of different types.
//' @param radius Radius of interaction.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
[[nodiscard]] Rcpp::NumericVector compute_delta_phi_dispersion(Rcpp::S4 configuration, Rcpp::NumericVector location, R_xlen_t type, int number_types, double radius = 0) {
  if(radius == 0) {
    return compute_delta_phi_dispersion_helper<Varphi::identity>(configuration, location, type, number_types);
  } else {
    return compute_delta_phi_dispersion_helper<Varphi::Strauss>(configuration, location, type, number_types, radius * radius);
  }
}

