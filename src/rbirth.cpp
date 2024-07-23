#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_wrapper.hpp"
#include "configuration/make_R_configuration.hpp"

#include "saturated_model/exponential_family_model.hpp"

#include "simulation/birth_death.hpp"
#include "simulation/rstrat_single.hpp"

#include "utility/construct_if_missing.hpp"
#include "utility/get_list_or_first_element.hpp"
#include "utility/get_number_types.hpp"
#include "utility/make_default_types.hpp"
#include "utility/window.hpp"

#include <random>
#include <vector>

// [[Rcpp::export]]
SEXP rbirth_cpp(R_xlen_t nsim,
                double horizon,
                Rcpp::CharacterVector types,
                R_xlen_t nquad,
                SEXP window,
                bool drop,
                R_xlen_t seed,
                Rcpp::NumericVector mark_range,
                SEXP starting_configuration,
                SEXP dummy,
                int nthreads,
                SEXP birth_alpha,
                Rcpp::NumericVector birth_beta0,
                SEXP birth_covariates,
                SEXP birth_beta,
                SEXP birth_gamma,
                SEXP birth_short_range,
                SEXP birth_medium_range,
                SEXP birth_long_range,
                R_xlen_t birth_saturation,
                SEXP birth_model,
                SEXP birth_medium_range_model,
                SEXP death_alpha,
                Rcpp::NumericVector death_beta0,
                SEXP death_covariates,
                SEXP death_beta,
                SEXP death_gamma,
                SEXP death_short_range,
                SEXP death_medium_range,
                SEXP death_long_range,
                R_xlen_t death_saturation,
                SEXP death_model,
                SEXP death_medium_range_model) {
  using Configuration_type = std::vector<ppjsdm::Marked_point>;

  // Compute number of types
  const auto number_types(types.size());

  // Convert starting configuration into C++ format
  Configuration_type vector_starting_configuration{};

  if(starting_configuration != R_NilValue) {
    const ppjsdm::Configuration_wrapper wrapped_configuration(Rcpp::wrap(starting_configuration));
    const auto length_configuration(ppjsdm::size(wrapped_configuration));
    vector_starting_configuration = Configuration_type(length_configuration);

    for(decltype(ppjsdm::size(wrapped_configuration)) j(0); j < length_configuration; ++j) {
      vector_starting_configuration[j] = wrapped_configuration[j];
    }
  }

  // Convert dummy into C++ format
  const ppjsdm::Configuration_wrapper wrapped_dummy(Rcpp::wrap(dummy));
  const auto length_dummy(ppjsdm::size(wrapped_dummy));
  auto vector_dummy = Configuration_type(length_dummy);

  for(decltype(ppjsdm::size(wrapped_dummy)) j(0); j < length_dummy; ++j) {
    vector_dummy[j] = wrapped_dummy[j];
  }

  std::mt19937 generator(seed);
  const auto cpp_window(ppjsdm::Window(window, mark_range));

  // Make the birth model
  const ppjsdm::Truncated_exponential_family_model_over_window<std::vector<double>> birth_exponential_model(cpp_window,
                                                                                                            Rcpp::as<std::vector<double>>(birth_beta0),
                                                                                                            birth_model,
                                                                                                            birth_medium_range_model,
                                                                                                            birth_alpha,
                                                                                                            birth_beta,
                                                                                                            birth_gamma,
                                                                                                            birth_covariates,
                                                                                                            birth_short_range,
                                                                                                            birth_medium_range,
                                                                                                            birth_long_range,
                                                                                                            birth_saturation);

  // Make the death model
  const ppjsdm::Truncated_exponential_family_model_over_window<std::vector<double>> death_exponential_model(cpp_window,
                                                                                                            Rcpp::as<std::vector<double>>(death_beta0),
                                                                                                            death_model,
                                                                                                            death_medium_range_model,
                                                                                                            death_alpha,
                                                                                                            death_beta,
                                                                                                            death_gamma,
                                                                                                            death_covariates,
                                                                                                            death_short_range,
                                                                                                            death_medium_range,
                                                                                                            death_long_range,
                                                                                                            death_saturation);

  Rcpp::List samples(nsim);
  Rcpp::List number_points(nsim);

  for(R_xlen_t i(0); i < nsim; ++i) {
    const auto sample(ppjsdm::generate_birth_death(generator,
                                                   cpp_window,
                                                   vector_dummy,
                                                   number_types,
                                                   horizon,
                                                   vector_starting_configuration,
                                                   birth_exponential_model,
                                                   death_exponential_model,
                                                   nthreads));
    samples[i] = ppjsdm::make_R_configuration(sample.first, types);
    Rcpp::as<Rcpp::List>(samples[i]).attr("number_points") = sample.second;
  }
  return ppjsdm::get_list_or_first_element(samples, nsim == 1 && drop);
}
