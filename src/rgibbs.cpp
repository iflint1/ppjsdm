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

template<typename Configuration, typename Model>
inline auto sample(const Model& model) {
  return ppjsdm::simulate_coupling_from_the_past<Configuration>(model);
}

template<typename Configuration, typename Model>
inline auto sample(const Model& model, R_xlen_t steps) {
  return ppjsdm::simulate_metropolis_hastings<Configuration>(model, steps);
}

template<typename Configuration, typename Model>
inline auto sample(const Model& model, R_xlen_t steps, const Configuration& starting_configuration) {
  return ppjsdm::simulate_metropolis_hastings<Configuration>(model, steps, starting_configuration);
}

template<typename Configuration, typename... Args>
inline SEXP rgibbs_helper(R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, const Args&... args) {
  Rcpp::List samples(nsim);

  for(R_xlen_t i(0); i < nsim; ++i) {
    const auto sample(sample<Configuration>(args...));
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
                Rcpp::NumericVector mark_range,
                SEXP starting_configuration) {
  using Configuration_type = std::vector<ppjsdm::Marked_point>;

  const auto cpp_window(ppjsdm::Window(window, mark_range));
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

  if(steps == 0) {
    return rgibbs_helper<Configuration_type>(nsim, types, drop, exponential_model);
  } else {
    if(starting_configuration != R_NilValue) {
      const ppjsdm::Configuration_wrapper wrapped_configuration(Rcpp::wrap(starting_configuration));
      const auto length_configuration(ppjsdm::size(wrapped_configuration));
      Configuration_type vector_starting_configuration(length_configuration);
      for(decltype(ppjsdm::size(wrapped_configuration)) j(0); j < length_configuration; ++j) {
        vector_starting_configuration[j] = wrapped_configuration[j];
      }
      return rgibbs_helper<Configuration_type>(nsim, types, drop, exponential_model, steps, vector_starting_configuration);
    } else {
      return rgibbs_helper<Configuration_type>(nsim, types, drop, exponential_model, steps);
    }
  }
}
