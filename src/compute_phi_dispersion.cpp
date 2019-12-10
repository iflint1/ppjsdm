#include <Rcpp.h>

#include "compute_phi_dispersion.h"

//' Compute the delta of the phi-dispersion of a marked configuration.
//'
//' @param configuration The configuration.
//' @param location Point to be added.
//' @param type Type of point to be added.
//' @param number_types Number of different types.
//' @param model String representing the model to simulate from. At the moment, either "i", "s" or "g".
//' @param radius Radius of interaction.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
[[nodiscard]] Rcpp::NumericVector compute_delta_phi_dispersion(Rcpp::S4 configuration, Rcpp::NumericVector location, R_xlen_t type, int number_types, Rcpp::CharacterVector model = "i", double radius = 0) {
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
  } else if(model[0] == "g"){
    const auto varphi{Geyer_papangelou{radius * radius, 2.0}};
    return varphi.compute(Rcpp::NumericVector(configuration.slot("x")),
                          Rcpp::NumericVector(configuration.slot("y")),
                          Rcpp::IntegerVector(configuration.slot("types")),
                          location, type, number_types);
  } else {
    Rcpp::Rcerr << "Incorrect model entered.\n";
    return {};
  }
}
