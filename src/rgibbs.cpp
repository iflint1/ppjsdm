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

#include <utility> // std::forward
#include <vector> // std::vector

namespace detail {

template<typename Configuration, typename... Args>
inline Configuration sample(Args&&... args);

template<typename Configuration, typename Model>
inline auto sample(const Model& model) {
  return ppjsdm::simulate_coupling_from_the_past<Configuration>(model);
}

template<typename Configuration, typename Model, typename... Args>
inline auto sample(const Model& model, R_xlen_t steps, Args&&... args) {
  return ppjsdm::simulate_metropolis_hastings<Configuration>(model, steps, std::forward<Args>(args)...);
}

} // namespace detail

template<typename Configuration, typename... Args>
inline SEXP rgibbs_helper(R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, const Args&... args) {
  Rcpp::List samples(nsim);

  for(R_xlen_t i(0); i < nsim; ++i) {
    const auto sample(detail::sample<Configuration>(args...));
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
                SEXP model,
                SEXP medium_range_model,
                bool drop,
                Rcpp::NumericVector mark_range,
                Rcpp::IntegerVector only_simulate_these_types,
                SEXP conditional_configuration,
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
    // TODO: Do not discard only_simulate_these_types or conditional_configuration.
    return rgibbs_helper<Configuration_type>(nsim, types, drop, exponential_model);
  } else {
    Configuration_type vector_conditional_configuration{};
    if(conditional_configuration != R_NilValue) {
      const ppjsdm::Configuration_wrapper wrapped_conditional_configuration(Rcpp::wrap(conditional_configuration));
      const auto length_conditional_configuration(ppjsdm::size(wrapped_conditional_configuration));
      vector_conditional_configuration = Configuration_type(length_conditional_configuration);
      for(decltype(ppjsdm::size(wrapped_conditional_configuration)) j(0); j < length_conditional_configuration; ++j) {
        vector_conditional_configuration[j] = wrapped_conditional_configuration[j];
      }
    }

    if(starting_configuration != R_NilValue) {
      const ppjsdm::Configuration_wrapper wrapped_configuration(Rcpp::wrap(starting_configuration));
      const auto length_configuration(ppjsdm::size(wrapped_configuration));
      Configuration_type vector_starting_configuration(length_configuration);
      for(decltype(ppjsdm::size(wrapped_configuration)) j(0); j < length_configuration; ++j) {
        vector_starting_configuration[j] = wrapped_configuration[j];
      }
      return rgibbs_helper<Configuration_type>(nsim, types, drop, exponential_model, steps, vector_starting_configuration, only_simulate_these_types, vector_conditional_configuration);
    } else {
      return rgibbs_helper<Configuration_type>(nsim, types, drop, exponential_model, steps, only_simulate_these_types, vector_conditional_configuration);
    }
  }
}
