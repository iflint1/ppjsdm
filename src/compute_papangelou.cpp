#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_wrapper.hpp"
#include "configuration/get_number_points.hpp"

#include "point/point_manipulation.hpp"

#include "saturated_model/exponential_family_model.hpp"

// [[Rcpp::export]]
Rcpp::NumericVector compute_papangelou_cpp(SEXP configuration,
                                           Rcpp::NumericVector x,
                                           Rcpp::NumericVector y,
                                           Rcpp::CharacterVector model,
                                           Rcpp::CharacterVector medium_range_model,
                                           Rcpp::List alpha,
                                           Rcpp::NumericVector beta0,
                                           Rcpp::NumericMatrix beta,
                                           Rcpp::NumericMatrix gamma,
                                           Rcpp::List covariates,
                                           Rcpp::List short_range,
                                           Rcpp::NumericMatrix medium_range,
                                           Rcpp::NumericMatrix long_range,
                                           R_xlen_t saturation,
                                           Rcpp::IntegerVector type,
                                           Rcpp::NumericVector mark,
                                           int nthreads) {
  const ppjsdm::Configuration_wrapper wrapped_configuration(configuration);
  std::vector<ppjsdm::Marked_point> points(x.size());
  for(typename decltype(points)::size_type i(0); i < points.size(); ++i) {
    points[i] = ppjsdm::Marked_point(x[i], y[i], type[i] - 1, mark[i]);
  }

  const ppjsdm::Truncated_exponential_family_model<Rcpp::NumericVector> exponential_model(beta0, model,
                                                                                          medium_range_model,
                                                                                          alpha, beta, gamma,
                                                                                          covariates, short_range,
                                                                                          medium_range, long_range,
                                                                                          saturation);
  return Rcpp::wrap(exponential_model.compute_papangelou_vectorized(points, wrapped_configuration, nthreads));
}
