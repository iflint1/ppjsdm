#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_wrapper.hpp"
#include "configuration/get_number_points.hpp"

#include "point/point_manipulation.hpp"

#include "saturated_model/exponential_family_model.hpp"

// [[Rcpp::export]]
Rcpp::NumericVector compute_papangelou_cpp(SEXP configuration, Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::CharacterVector model, Rcpp::CharacterVector medium_range_model, Rcpp::NumericMatrix alpha, Rcpp::NumericVector beta0, Rcpp::NumericMatrix beta, Rcpp::NumericMatrix gamma, Rcpp::List covariates, Rcpp::NumericMatrix short_range, Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range, R_xlen_t saturation, R_xlen_t type = 1, double mark = 1.0) {
  const ppjsdm::Configuration_wrapper wrapped_configuration(configuration);
  const ppjsdm::Truncated_exponential_family_model<Rcpp::NumericVector> exponential_model(beta0, model, medium_range_model, alpha, beta, gamma, covariates, short_range, medium_range, long_range, saturation);
  const auto length_x(x.size());
  std::vector<ppjsdm::Marked_point> points(x.size());
  for(typename decltype(points)::size_type i(0); i < points.size(); ++i) {
    points[i] = ppjsdm::Marked_point(x[i], y[i], type - 1, mark);
  }
  const auto result_vector(exponential_model.compute_papangelou_vectorized(points, wrapped_configuration));
  Rcpp::NumericVector result(Rcpp::no_init(length_x));
  for(R_xlen_t i(0); i < length_x; ++i) {
    result[i] = result_vector[i];
  }
  return result;
}
