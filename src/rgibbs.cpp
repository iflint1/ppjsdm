#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/make_R_configuration.hpp"

#include "point/point_manipulation.hpp"

#include "saturated_varphi_model/model.hpp"

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
SEXP rgibbs_cpp(SEXP window, SEXP alpha, SEXP lambda, SEXP covariates, SEXP coefs, SEXP radius, R_xlen_t saturation, R_xlen_t steps, R_xlen_t nsim, SEXP types, Rcpp::CharacterVector model, bool drop) {
  if(Rf_isNull(covariates)) {
    covariates = Rcpp::wrap(Rcpp::List(0));
  }

  const auto number_types(ppjsdm::get_number_types_and_check_conformance(alpha, lambda, radius, types));
  alpha = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(alpha, 0., number_types);
  lambda = ppjsdm::construct_if_missing<Rcpp::NumericVector>(lambda, 1., number_types);

  const auto coefs_nrows(number_types);
  const auto coefs_ncols(Rcpp::as<Rcpp::List>(covariates).size());
  coefs = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(coefs, 1., coefs_nrows, coefs_ncols);
  if(Rcpp::as<Rcpp::NumericMatrix>(coefs).nrow() != coefs_nrows || Rcpp::as<Rcpp::NumericMatrix>(coefs).ncol() != coefs_ncols) {
    Rcpp::stop("The parameter `coefs` does not have the right dimensions.");
  }

  types = ppjsdm::make_types(types, number_types, lambda);
  // TODO: Think about what format for coefs, see also compute_papangelou.cpp.
  return ppjsdm::call_on_wrapped_window(window, [alpha, lambda, coefs, covariates, radius, saturation, steps, nsim, types, model, drop, number_types](const auto& w) {
    return ppjsdm::call_on_list_or_vector(lambda, [alpha, lambda, coefs, covariates, radius, saturation, steps, nsim, types, model, drop, number_types, &w](const auto& l) {
      const auto r(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(radius, 0.1 * w.diameter(), number_types));
      return ppjsdm::call_on_model(model, alpha, l, coefs, covariates, r, saturation, [&w, steps, nsim, types, drop, number_types](const auto& model) {
        if(steps == 0) {
          return rgibbs_helper(model, w, nsim, types, drop, number_types);
        } else {
          return rgibbs_helper(model, w, nsim, types, drop, number_types, steps);
        }
      });
    });
  });
}
