#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/make_R_configuration.hpp"

#include "phi_dispersion_model/compute_phi_dispersion.hpp"

#include "point/point_manipulation.hpp"

#include "simulation/coupling_from_the_past.hpp"
#include "simulation/metropolis_hastings.hpp"

#include "utility/call_on_list_or_vector.hpp"
#include "utility/construct_if_missing.hpp"
#include "utility/get_list_or_first_element.hpp"
#include "utility/get_number_types.hpp"
#include "utility/make_default_types.hpp"
#include "utility/window_utilities.hpp"

#include <vector> // std::vector

template<typename Model, typename Window>
inline SEXP rgibbs_helper(const Model& model, const Window& window, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, R_xlen_t point_types) {
  Rcpp::List samples(nsim);

  for(R_xlen_t i(0); i < nsim; ++i) {
    const auto sample(ppjsdm::simulate_coupling_from_the_past<std::vector<ppjsdm::Marked_point>>(model, window, point_types));
    samples[i] = ppjsdm::make_R_configuration(sample, types);
  }

  return ppjsdm::get_list_or_first_element(samples, nsim, drop);
}

template<typename Model, typename Window>
inline SEXP rgibbs_helper(const Model& model, const Window& window, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, R_xlen_t point_types, R_xlen_t steps) {
  Rcpp::List samples(nsim);

  for(R_xlen_t i(0); i < nsim; ++i) {
    const auto sample(ppjsdm::simulate_metropolis_hastings<std::vector<ppjsdm::Marked_point>>(model, window, steps, point_types));
    samples[i] = ppjsdm::make_R_configuration(sample, types);
  }

  return ppjsdm::get_list_or_first_element(samples, nsim, drop);
}

// [[Rcpp::export]]
SEXP rgibbs_cpp(SEXP window, SEXP alpha, SEXP lambda, SEXP covariates, SEXP coefs, SEXP radius, R_xlen_t steps, R_xlen_t nsim, SEXP types, Rcpp::CharacterVector model, bool drop) {
  if(Rf_isNull(covariates)) {
    covariates = Rcpp::wrap(Rcpp::List(0));
  }

  const auto point_types(ppjsdm::get_number_types_and_check_conformance(alpha, lambda, radius, types));
  alpha = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(point_types, alpha, 0.);
  lambda = ppjsdm::construct_if_missing<Rcpp::NumericVector>(point_types, lambda, 1.);
  radius = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(point_types, radius, 0.);
  coefs = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(Rcpp::as<Rcpp::List>(covariates).size(), coefs, 0.);
  types = ppjsdm::make_types(types, point_types, lambda);
  // TODO: Think about what format for coefs.
  return ppjsdm::call_on_wrapped_window(window, [alpha, lambda, coefs, covariates, radius, steps, nsim, types, model, drop, point_types](const auto& w) {
    return ppjsdm::call_on_list_or_vector(lambda, [alpha, lambda, coefs, covariates, radius, steps, nsim, types, model, drop, point_types, &w](const auto& l) {
      return ppjsdm::call_on_model(model, alpha, l, coefs, covariates, radius, [&w, steps, nsim, types, drop, point_types](const auto& model){
        if(steps == 0) {
          return rgibbs_helper(model, w, nsim, types, drop, point_types);
        } else {
          return rgibbs_helper(model, w, nsim, types, drop, point_types, steps);
        }
      });
    });
  });
}
