#include <Rcpp.h>

#include <vector> // std::vector

#include "compute_phi_dispersion.h"
#include "configuration_utilities.h"

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
  const ppjsdm::Configuration_wrapper wrapped_configuration(configuration);
  const auto number_points(wrapped_configuration.get_number_points());
  if(model[0] == "identity") {
    const ppjsdm::Varphi_model_papangelou<ppjsdm::varphi::Identity> varphi{};
    return varphi.compute(wrapped_configuration.x(), wrapped_configuration.y(), wrapped_configuration.types(), location, type, number_types, number_points);
  } else if(model[0] == "Strauss") {
    const ppjsdm::Varphi_model_papangelou<ppjsdm::varphi::Strauss> varphi(radius);
    return varphi.compute(wrapped_configuration.x(), wrapped_configuration.y(), wrapped_configuration.types(), location, type, number_types, number_points);
  } else if(model[0] == "Geyer"){
    const ppjsdm::Geyer_papangelou varphi(radius, 2.0);
    return varphi.compute(wrapped_configuration.x(), wrapped_configuration.y(), wrapped_configuration.types(), location, type, number_types, number_points);
  } else if(model[0] == "neighbour"){
    const ppjsdm::Nearest_neighbour_papangelou<ppjsdm::varphi::Identity> varphi{};
    // TODO: At the moment, I'm using push_backs in the computation, so I need to convert to std::vector first
    std::vector<double> x_vector(number_points);
    std::vector<double> y_vector(number_points);
    std::vector<int> types_vector(number_points);
    for(R_xlen_t i(0); i < number_points; ++i) {
      x_vector[i] = wrapped_configuration.x(i);
      y_vector[i] = wrapped_configuration.y(i);
      types_vector[i] = wrapped_configuration.types(i);
    }
    return varphi.compute(x_vector,
                          y_vector,
                          types_vector,
                          location, type, number_types, number_points);
  } else {
    Rcpp::stop("Incorrect model entered.\n");
  }
}
