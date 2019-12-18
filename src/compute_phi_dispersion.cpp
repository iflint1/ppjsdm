#include <Rcpp.h>

#include "compute_phi_dispersion.h"

//' Compute the delta of the phi-dispersion of a marked configuration.
//'
//' @param configuration The configuration.
//' @param location Point to be added.
//' @param type Type of point to be added.
//' @param number_types Number of different types.
//' @param model String representing the model to simulate from. At the moment, either "identity", "Strauss" or "Geyer".
//' @param radius Radius of interaction.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
[[nodiscard]] Rcpp::NumericVector compute_delta_phi_dispersion(Rcpp::S4 configuration, Rcpp::NumericVector location, R_xlen_t type, int number_types, Rcpp::CharacterVector model = "identity", double radius = 0) {
  const auto x{Rcpp::NumericVector(configuration.slot("x"))};
  if(model[0] == "identity") {
    const auto varphi{Varphi_model_papangelou<varphi::Identity>{}};
    return varphi.compute(x,
                          Rcpp::NumericVector(configuration.slot("y")),
                          Rcpp::IntegerVector(configuration.slot("types")),
                          location, type, number_types, x.size());
  } else if(model[0] == "Strauss") {
    const auto varphi{Varphi_model_papangelou<varphi::Strauss>{radius}};
    return varphi.compute(x,
                          Rcpp::NumericVector(configuration.slot("y")),
                          Rcpp::IntegerVector(configuration.slot("types")),
                          location, type, number_types, x.size());
  } else if(model[0] == "Geyer"){
    const auto varphi{Geyer_papangelou{radius, 2.0}};
    return varphi.compute(x,
                          Rcpp::NumericVector(configuration.slot("y")),
                          Rcpp::IntegerVector(configuration.slot("types")),
                          location, type, number_types, x.size());
  } else {
    Rcpp::stop("Incorrect model entered.\n");
    return {};
  }
}
