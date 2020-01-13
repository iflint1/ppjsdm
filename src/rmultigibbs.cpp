#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/make_R_configuration.h"

#include "phi_dispersion_model/compute_phi_dispersion.h"

#include "point/point_manipulation.h"

#include "simulation/simulate_metropolis_hastings.h"

#include "utility/call_on_list_or_vector.h"
#include "utility/construct_if_missing.h"
#include "utility/get_list_or_first_element.h"
#include "utility/get_number_types.h"
#include "utility/make_default_types.h"
#include "utility/window_utilities.h"

#include <vector> // std::vector

template<typename Model, typename Window>
inline SEXP rmultigibbs_helper(const Model& model, const Window& window, R_xlen_t steps, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, R_xlen_t point_types) {
  Rcpp::List samples(nsim);

  for(R_xlen_t i(0); i < nsim; ++i) {
    const auto sample(ppjsdm::simulate_metropolis_hastings<std::vector<ppjsdm::Marked_point>>(model, window, steps, point_types));
    samples[i] = ppjsdm::make_R_configuration(sample, types);
  }

  return ppjsdm::get_list_or_first_element(samples, nsim, drop);
}

//' Sample a multivariate Gibbs point processes
//'
//' @param window Simulation window.
//' @param alpha Repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.
//' @param lambda A vector representing the intensities of the point processes.
//' Default is a vector of same size as types, filled with ones.
//' @param nu A vector representing the dispersion of the number of points.
//' Default is a vector of same size as types, filled with ones.
//' @param radius Symmetric matrix of interaction radii. Filled by zeroes by default;
//' @param steps Number of steps in the Metropolis algorithm. Default is 30000.
//' @param nsim Number of samples to generate. Default is 1.
//' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
//' @param model String representing the model to simulate from. At the moment, either "identity", "Strauss", "Geyer" or "neighbour",
//' with default beign "identity".
//' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration. Default is TRUE.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
SEXP rmultigibbs(SEXP window = R_NilValue, SEXP alpha = R_NilValue, SEXP lambda = R_NilValue, SEXP nu = R_NilValue, SEXP radius = R_NilValue, R_xlen_t steps = 30000, R_xlen_t nsim = 1, SEXP types = R_NilValue, Rcpp::CharacterVector model = "identity", bool drop = true) {
  const auto point_types(ppjsdm::get_number_types_and_check_conformance(alpha, lambda, nu, radius, types));
  alpha = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(point_types, alpha, 0.);
  lambda = ppjsdm::construct_if_missing<Rcpp::NumericVector>(point_types, lambda, 1.);
  nu = ppjsdm::construct_if_missing<Rcpp::NumericVector>(point_types, nu, 1.);
  radius = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(point_types, radius, 0.);
  types = ppjsdm::make_types(types, point_types, lambda, nu);
  return ppjsdm::call_on_wrapped_window(window, [&alpha, &lambda, &nu, &radius, steps, nsim, &types, &model, drop, point_types](const auto& w) {
    return ppjsdm::call_on_list_or_vector(lambda, [&alpha, &lambda, &nu, &radius, steps, nsim, &types, &model, drop, point_types, &w](const auto& l) {
      return ppjsdm::call_on_list_or_vector(nu, [&alpha, &lambda, &nu, &radius, steps, nsim, &types, &model, drop, point_types, &w, &l](const auto& n) {
        return ppjsdm::call_on_model(model, alpha, l, n, radius, [&w, steps, nsim, &types, drop, point_types](const auto& model){
          return rmultigibbs_helper(model, w, steps, nsim, types, drop, point_types);
        });
      });
    });
  });
}
