#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_wrapper.hpp"
#include "configuration/get_number_points.hpp"

#include "point/point_manipulation.hpp"

#include "saturated_varphi_model/model.hpp"

//' Compute Papangelou conditional intensity of the model.
//'
//' @param configuration Configuration.
//' @param coordinates Coordinates of the point at which to evaluate the Papangelou conditional intensity.
//' @param type Type of the point (as an integer >= 1).
//' @param model String representing the model to simulate from. You can check the currently authorised models with a call to `show_model()`.
//' @param alpha Repulsion matrix.
//' @param lambda A vector representing the intensities of the point processes.
//' @param coefs Coefficients related to covariates.
//' @param covariates Covariates.
//' @param radius Interaction radii.
//' @param saturation Saturation parameter.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
double compute_papangelou(SEXP configuration, Rcpp::NumericVector coordinates, R_xlen_t type, Rcpp::CharacterVector model, Rcpp::NumericMatrix alpha, Rcpp::NumericVector lambda, Rcpp::NumericMatrix coefs, Rcpp::List covariates, Rcpp::NumericMatrix radius, R_xlen_t saturation) {
  const ppjsdm::Configuration_wrapper wrapped_configuration(configuration);
  return ppjsdm::call_on_model(model, alpha, lambda, coefs, covariates, radius, saturation, [&wrapped_configuration, &coordinates, type](const auto& model) {
    const auto points_by_type(ppjsdm::get_number_points(wrapped_configuration));
    const auto number_types(points_by_type.size());
    return model.compute_papangelou(wrapped_configuration, ppjsdm::Marked_point(coordinates[0], coordinates[1], type - 1), number_types);
  });
}
