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
// [[nodiscard]] Rcpp::NumericMatrix compute_phi_dispersion(Rcpp::S4 configuration, int number_types, double radius = 0) {
//   if(radius == 0) {
//     const auto varphi{varphi::Identity{}};
//     return compute_phi_dispersion_helper(configuration, number_types, varphi);
//   } else {
//     const auto varphi{varphi::Strauss{radius * radius}};
//     return compute_phi_dispersion_helper(configuration, number_types, varphi);
//   }
// }

//' Compute the delta of the phi-dispersion of a marked configuration (radius = 0 corresponds to the proposal model).
//'
//' @param configuration The configuration.
//' @param location Point to be added.
//' @param type Type of point to be added.
//' @param number_types Number of different types.
//' @param model String representing the model to simulate from.
//' @param radius Radius of interaction.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
[[nodiscard]] Rcpp::NumericVector compute_delta_phi_dispersion(Rcpp::S4 configuration, Rcpp::NumericVector location, R_xlen_t type, int number_types, Rcpp::CharacterVector model, double radius = 0) {
  // TODO: Error message on incorrect model specification.
  if(model[0] == "i") {
    const auto varphi{Varphi_model_papangelou<varphi::Identity>{}};
    return varphi.compute(Rcpp::NumericVector(configuration.slot("x")),
                          Rcpp::NumericVector(configuration.slot("y")),
                          Rcpp::IntegerVector(configuration.slot("types")),
                          location, type, number_types);
  } else if(model[0] == "s") {
    const auto varphi{Varphi_model_papangelou<varphi::Strauss>{radius * radius}};
    return varphi.compute(Rcpp::NumericVector(configuration.slot("x")),
                          Rcpp::NumericVector(configuration.slot("y")),
                          Rcpp::IntegerVector(configuration.slot("types")),
                          location, type, number_types);
  } else {
    const auto varphi{Geyer_papangelou{radius * radius, 2.0}};
    return varphi.compute(Rcpp::NumericVector(configuration.slot("x")),
                          Rcpp::NumericVector(configuration.slot("y")),
                          Rcpp::IntegerVector(configuration.slot("types")),
                          location, type, number_types);
  }
}

//' Compute the delta of the phi-dispersion of a marked configuration in the Geyer model setting.
//'
//' @param configuration The configuration.
//' @param location Point to be added.
//' @param type Type of point to be added.
//' @param number_types Number of different types.
//' @param radius Radius of interaction.
//' @param saturation Saturation parameter.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
//[[nodiscard]] Rcpp::NumericVector compute_delta_phi_dispersion_geyer(Rcpp::S4 configuration, Rcpp::NumericVector location, R_xlen_t type, int number_types, double radius, double saturation) {
//  return compute_delta_phi_dispersion_geyer_helper(configuration, location, type, number_types, radius * radius, saturation);
//}

