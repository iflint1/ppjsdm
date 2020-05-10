#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/make_R_configuration.hpp"

#include "point/point_manipulation.hpp"

#include "saturated_model/exponential_family_model.hpp"

#include "simulation/coupling_from_the_past.hpp"
#include "simulation/metropolis_hastings.hpp"

#include "utility/construct_if_missing.hpp"
#include "utility/get_list_or_first_element.hpp"
#include "utility/get_number_types.hpp"
#include "utility/is_symmetric_matrix.hpp"
#include "utility/make_default_types.hpp"
#include "utility/window.hpp"

#include <vector> // std::vector

template<typename Model>
inline SEXP rgibbs_helper(const Model& model, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, R_xlen_t number_types) {
  Rcpp::List samples(nsim);

  for(R_xlen_t i(0); i < nsim; ++i) {
    const auto sample(ppjsdm::simulate_coupling_from_the_past<std::vector<ppjsdm::Marked_point>>(model, number_types));
    samples[i] = ppjsdm::make_R_configuration(sample, types);
  }

  return ppjsdm::get_list_or_first_element(samples, nsim == 1 && drop);
}

template<typename Model>
inline SEXP rgibbs_helper(const Model& model, const ppjsdm::Window& window, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, R_xlen_t number_types, R_xlen_t steps) {
  Rcpp::List samples(nsim);

  for(R_xlen_t i(0); i < nsim; ++i) {
    const auto sample(ppjsdm::simulate_metropolis_hastings<std::vector<ppjsdm::Marked_point>>(model, window, steps, number_types));
    samples[i] = ppjsdm::make_R_configuration(sample, types);
  }

  return ppjsdm::get_list_or_first_element(samples, nsim == 1 && drop);
}

// [[Rcpp::export]]
SEXP rgibbs_cpp(SEXP window, SEXP alpha, Rcpp::NumericVector beta0, SEXP covariates, SEXP beta, SEXP gamma, SEXP short_range, SEXP medium_range, SEXP long_range, R_xlen_t saturation, R_xlen_t steps, R_xlen_t nsim, SEXP types, Rcpp::CharacterVector model, Rcpp::CharacterVector medium_range_model, bool drop, Rcpp::NumericVector mark_range) {
  if(Rf_isNull(covariates)) {
    covariates = Rcpp::wrap(Rcpp::List(0));
  }

  // TODO: Centralize defaults
  const auto number_types(ppjsdm::get_number_types_and_check_conformance(alpha, gamma, beta0, short_range, medium_range, long_range, types));
  alpha = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(alpha, 0., number_types);
  gamma = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(gamma, 0., number_types);
  beta0 = ppjsdm::construct_if_missing<Rcpp::NumericVector>(beta0, 0., number_types);
  if(!ppjsdm::is_symmetric_matrix(alpha) || !ppjsdm::is_symmetric_matrix(gamma)) {
    Rcpp::stop("Either alpha or gamma is not symmetric.");
  }

  const auto beta_nrows(number_types);
  const auto beta_ncols(Rcpp::as<Rcpp::List>(covariates).size());
  beta = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(beta, 1., beta_nrows, beta_ncols);
  if(beta_ncols != 0 && (Rcpp::as<Rcpp::NumericMatrix>(beta).nrow() != beta_nrows
                            || Rcpp::as<Rcpp::NumericMatrix>(beta).ncol() != beta_ncols)) {
    Rcpp::stop("The parameter `beta` does not have the right dimensions.");
  }

  types = ppjsdm::make_types(types, number_types, beta0);
  const auto cpp_window(ppjsdm::Window(window, mark_range));
  const auto sh(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(short_range, 0.1 * cpp_window.diameter(), number_types));
  const auto me(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(medium_range, 0., number_types));
  const auto lo(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(long_range, 0., number_types));
  if(!ppjsdm::is_symmetric_matrix(sh) || !ppjsdm::is_symmetric_matrix(me) || !ppjsdm::is_symmetric_matrix(lo)) {
    Rcpp::stop("One of the interaction radii matrices is not symmetric.");
  }
  if(steps == 0) {
    const ppjsdm::Truncated_exponential_family_model_over_window<Rcpp::NumericVector> exponential_model(cpp_window, beta0, model, medium_range_model, alpha, beta, gamma, covariates, sh, me, lo, saturation);
    return rgibbs_helper(exponential_model, nsim, types, drop, number_types);

  } else {
    const ppjsdm::Truncated_exponential_family_model<Rcpp::NumericVector> exponential_model(beta0, model, medium_range_model, alpha, beta, gamma, covariates, sh, me, lo, saturation);
    return rgibbs_helper(exponential_model, cpp_window, nsim, types, drop, number_types, steps);
  }
}
