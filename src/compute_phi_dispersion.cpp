#include <Rcpp.h>

#include <vector> // std::vector

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
Rcpp::NumericVector compute_delta_phi_dispersion(Rcpp::List configuration, Rcpp::NumericVector location, R_xlen_t type, int number_types, Rcpp::CharacterVector model = "identity", double radius = 0) {
  const auto x(Rcpp::as<Rcpp::NumericVector>(configuration["x"]));
  const auto y(Rcpp::as<Rcpp::NumericVector>(configuration["y"]));
  const auto types(Rcpp::as<Rcpp::IntegerVector>(configuration["types"]));
  if(model[0] == "identity") {
    const Varphi_model_papangelou<varphi::Identity> varphi{};
    return varphi.compute(x, y, types,
                          location, type, number_types, x.size());
  } else if(model[0] == "Strauss") {
    const Varphi_model_papangelou<varphi::Strauss> varphi(radius);
    return varphi.compute(x, y, types,
                          location, type, number_types, x.size());
  } else if(model[0] == "Geyer"){
    const Geyer_papangelou varphi(radius, 2.0);
    return varphi.compute(x, y, types,
                          location, type, number_types, x.size());
  } else if(model[0] == "neighbour"){
    const Nearest_neighbour_papangelou<varphi::Identity> varphi{};
    // TODO: At the moment, I'm using push_backs in the computation, so I need to convert to std::vector first
    std::vector<double> x_vector(x.size());
    std::vector<double> y_vector(x.size());
    std::vector<int> types_vector(x.size());
    for(int i{0}; i < x.size(); ++i) {
      x_vector[i] = x[i];
      y_vector[i] = y[i];
      types_vector[i] = types[i];
    }
    return varphi.compute(x_vector,
                          y_vector,
                          types_vector,
                          location, type, number_types, x.size());
  } else {
    Rcpp::stop("Incorrect model entered.\n");
  }
}
