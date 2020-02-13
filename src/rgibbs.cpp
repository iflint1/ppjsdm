#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/make_R_configuration.hpp"

#include "point/point_manipulation.hpp"

#include "saturated_varphi_model/exponential_family_model.hpp"

#include "simulation/coupling_from_the_past.hpp"
#include "simulation/metropolis_hastings.hpp"

#include "utility/call_on_list_or_vector.hpp"
#include "utility/construct_if_missing.hpp"
#include "utility/get_list_or_first_element.hpp"
#include "utility/get_number_types.hpp"
#include "utility/make_default_types.hpp"
#include "utility/window_utilities.hpp"

#include <vector> // std::vector

template<typename Model>
inline SEXP rgibbs_helper(const Model& model, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, R_xlen_t number_types) {
  Rcpp::List samples(nsim);

  for(R_xlen_t i(0); i < nsim; ++i) {
    const auto sample(ppjsdm::simulate_coupling_from_the_past<std::vector<ppjsdm::Marked_point>>(model, number_types));
    samples[i] = ppjsdm::make_R_configuration(sample, types);
  }

  return ppjsdm::get_list_or_first_element(samples, nsim, drop);
}

template<typename Model, typename Window>
inline SEXP rgibbs_helper(const Model& model, const Window& window, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, R_xlen_t number_types, R_xlen_t steps) {
  Rcpp::List samples(nsim);

  for(R_xlen_t i(0); i < nsim; ++i) {
    const auto sample(ppjsdm::simulate_metropolis_hastings<std::vector<ppjsdm::Marked_point>>(model, window, steps, number_types));
    samples[i] = ppjsdm::make_R_configuration(sample, types);
  }

  return ppjsdm::get_list_or_first_element(samples, nsim, drop);
}

// TODO: Make rgibbs work with R::Inf

// [[Rcpp::export]]
SEXP rgibbs_cpp(SEXP window, SEXP alpha, SEXP lambda, SEXP covariates, SEXP beta, SEXP gamma, SEXP short_range, SEXP medium_range, SEXP long_range, R_xlen_t saturation, R_xlen_t steps, R_xlen_t nsim, SEXP types, Rcpp::CharacterVector model, Rcpp::CharacterVector medium_range_model, bool drop) {
  if(Rf_isNull(covariates)) {
    covariates = Rcpp::wrap(Rcpp::List(0));
  }

  const auto number_types(ppjsdm::get_number_types_and_check_conformance(alpha, gamma, lambda, short_range, medium_range, long_range, types));
  alpha = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(alpha, 0., number_types);
  gamma = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(gamma, 0., number_types);
  lambda = ppjsdm::construct_if_missing<Rcpp::NumericVector>(lambda, 1., number_types);

  const auto beta_nrows(number_types);
  const auto beta_ncols(Rcpp::as<Rcpp::List>(covariates).size());
  beta = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(beta, 1., beta_nrows, beta_ncols);
  if(beta_ncols != 0 && (Rcpp::as<Rcpp::NumericMatrix>(beta).nrow() != beta_nrows
                            || Rcpp::as<Rcpp::NumericMatrix>(beta).ncol() != beta_ncols)) {
    Rcpp::stop("The parameter `beta` does not have the right dimensions.");
  }

  types = ppjsdm::make_types(types, number_types, lambda);
  return ppjsdm::call_on_wrapped_window(window, [alpha, lambda, beta, gamma, covariates, short_range, medium_range, long_range, saturation, steps, nsim, types, model, medium_range_model, drop, number_types](const auto& w) {
    return ppjsdm::call_on_list_or_vector(lambda, [alpha, lambda, beta, gamma, covariates, short_range, medium_range, long_range, saturation, steps, nsim, types, model, medium_range_model, drop, number_types, &w](const auto& l) {
      const auto sh(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(short_range, 0.1 * w.diameter(), number_types));
      const auto me(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(medium_range, 0.1 * w.diameter(), number_types));
      const auto lo(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(long_range, 0.2 * w.diameter(), number_types));
      if(steps == 0) {
        return ppjsdm::call_on_model(w, model, medium_range_model, l, sh, me, lo, saturation, [nsim, types, drop, number_types](const auto& model) {
          return rgibbs_helper(model, nsim, types, drop, number_types);
        }, alpha, beta, gamma, covariates);
      } else {
        return ppjsdm::call_on_model(model, medium_range_model, l, sh, me, lo, saturation, [&w, steps, nsim, types, drop, number_types](const auto& model) {
          return rgibbs_helper(model, w, nsim, types, drop, number_types, steps);
        }, alpha, beta, gamma, covariates);
      }
    });
  });
}
