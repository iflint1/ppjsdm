#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_wrapper.hpp"
#include "configuration/get_number_points.hpp"

#include "point/point_manipulation.hpp"

#include "saturated_model/exponential_family_model.hpp"

//' Compute Papangelou conditional intensity of the model.
//'
//' @param configuration Configuration.
//' @param x Coordinates along the x-axis of the points at which to evaluate the Papangelou conditional intensity.
//' @param y Coordinates along the x-axis of the points at which to evaluate the Papangelou conditional intensity.
//' @param type Type of the point (as an integer >= 1).
//' @param model String representing the model to use. You can check the currently authorised models with a call to `show_models()`.
//' @param medium_range_model String representing the medium range model to use. You can check the currently authorised models with a call to `show_medium_range_models()`.
//' @param alpha Short range repulsion matrix.
//' @param lambda A vector representing the intensities of the point processes.
//' @param beta Coefficients related to covariates.
//' @param gamma Medium range repulsion matrix.
//' @param covariates Covariates.
//' @param short_range Short range interaction radii.
//' @param medium_range Medium range interaction radii.
//' @param long_range Long range interaction radii.
//' @param saturation Saturation parameter.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
Rcpp::NumericVector compute_papangelou(SEXP configuration, Rcpp::NumericVector x, Rcpp::NumericVector y, R_xlen_t type, Rcpp::CharacterVector model, Rcpp::CharacterVector medium_range_model, Rcpp::NumericMatrix alpha, Rcpp::NumericVector lambda, Rcpp::NumericMatrix beta, Rcpp::NumericMatrix gamma, Rcpp::List covariates, Rcpp::NumericMatrix short_range, Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range, R_xlen_t saturation) {
  const ppjsdm::Configuration_wrapper wrapped_configuration(configuration);
  return ppjsdm::call_on_model(model, medium_range_model, lambda, short_range, medium_range, long_range, saturation, [&wrapped_configuration, &x, &y, type](const auto& model) {
    const auto length_x(x.size());
    Rcpp::NumericVector result(Rcpp::no_init(length_x));
    for(R_xlen_t i(0); i < length_x; ++i) {
      // TODO: Add mark to the parameters.
      result[i] = model.compute_papangelou(ppjsdm::Marked_point(x[i], y[i], type - 1, 1.0), wrapped_configuration);
    }
    return result;
  }, alpha, beta, gamma, covariates);
}
