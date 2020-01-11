#include <Rcpp.h>

#include "utility/compute_phi_dispersion.h"
#include "utility/configuration_manipulation.h"
#include "utility/configuration_wrappers.h"

//' Compute the delta of the phi-dispersion of a marked configuration.
//'
//' @param configuration Configuration.
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
  return ppjsdm::call_on_papangelou(model, radius, [&configuration, &location, type, number_types](const auto& papangelou) {
    const ppjsdm::Configuration_wrapper wrapped_configuration(configuration);
    const auto number_points(ppjsdm::size(wrapped_configuration));
    return papangelou.compute(wrapped_configuration, ppjsdm::Marked_point(location[0], location[1], type), number_types, number_points);
  });
}
