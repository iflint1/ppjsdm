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
SEXP rgibbs_cpp(SEXP window,
                SEXP alpha,
                Rcpp::NumericVector beta0,
                SEXP covariates,
                SEXP beta,
                SEXP gamma,
                SEXP short_range,
                SEXP medium_range,
                SEXP long_range,
                R_xlen_t saturation,
                R_xlen_t steps,
                R_xlen_t nsim,
                SEXP types,
                Rcpp::CharacterVector model,
                Rcpp::CharacterVector medium_range_model,
                bool drop,
                Rcpp::NumericVector mark_range) {
  const auto cpp_window(ppjsdm::Window(window, mark_range));
  const auto number_types(beta0.size());
  if(steps == 0) {
    const ppjsdm::Truncated_exponential_family_model_over_window<Rcpp::NumericVector> exponential_model(cpp_window,
                                                                                                        beta0,
                                                                                                        model,
                                                                                                        medium_range_model,
                                                                                                        alpha,
                                                                                                        beta,
                                                                                                        gamma,
                                                                                                        covariates,
                                                                                                        short_range,
                                                                                                        medium_range,
                                                                                                        long_range,
                                                                                                        saturation);
    return rgibbs_helper(exponential_model, nsim, types, drop, number_types);

  } else {
    const ppjsdm::Truncated_exponential_family_model<Rcpp::NumericVector> exponential_model(beta0,
                                                                                            model,
                                                                                            medium_range_model,
                                                                                            alpha,
                                                                                            beta,
                                                                                            gamma,
                                                                                            covariates,
                                                                                            short_range,
                                                                                            medium_range,
                                                                                            long_range,
                                                                                            saturation);
    return rgibbs_helper(exponential_model, cpp_window, nsim, types, drop, number_types, steps);
  }
}
