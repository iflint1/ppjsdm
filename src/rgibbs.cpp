#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/get_number_points.hpp"
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

#include <random> // Generator
#include <utility> // std::forward
#include <tuple> // std::pair
#include <vector> // std::vector

#ifdef _OPENMP
#include <omp.h>
#endif

namespace detail {

template<typename Configuration, typename Generator, typename... Args>
inline std::pair<Configuration, std::vector<decltype(get_number_points(Configuration{}))>> sample(Generator& generator, Args&&... args);

template<typename Configuration, typename Generator, typename Model>
inline auto sample(Generator& generator, const Model& model) {
  return ppjsdm::simulate_coupling_from_the_past<Configuration>(generator, model);
}

template<typename Configuration, typename Generator, typename Model, typename... Args>
inline auto sample(Generator& generator, const Model& model, R_xlen_t steps, bool debug, Args&&... args) {
  return ppjsdm::simulate_metropolis_hastings<Configuration>(generator, model, steps, debug, std::forward<Args>(args)...);
}

} // namespace detail

template<typename Configuration, typename... Args>
inline SEXP rgibbs_helper(R_xlen_t nthreads,
                          R_xlen_t seed,
                          R_xlen_t nsim,
                          Rcpp::CharacterVector types,
                          bool drop,
                          const Args&... args) {
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif

  Rcpp::List samples(nsim);
  std::vector<Configuration> cpp_samples(nsim);
  std::vector<std::vector<decltype(ppjsdm::get_number_points(Configuration{}))>> number_points(nsim);

#pragma omp parallel
{
  decltype(cpp_samples) cpp_samples_private(cpp_samples.size());
  decltype(number_points) number_points_private(number_points.size());
  R_xlen_t thread_num(0);
#ifdef _OPENMP
  thread_num = omp_get_thread_num();
#endif
  std::mt19937 generator(thread_num + seed);
#pragma omp for nowait
  for(typename decltype(cpp_samples_private)::size_type index = 0; index < cpp_samples_private.size(); ++index) {
    const auto sample(detail::sample<Configuration>(generator, args...));
    cpp_samples_private[index] = sample.first;
    number_points_private[index] = sample.second;
  }
#pragma omp critical
  for(std::remove_cv_t<decltype(cpp_samples_private.size())> index(0); index < cpp_samples_private.size(); ++index) {
    if(cpp_samples_private[index] != Configuration{} ||
       number_points_private[index] != typename decltype(number_points)::value_type{}) {
      cpp_samples[index] = cpp_samples_private[index];
      number_points[index] = number_points_private[index];
    }
  }
}

  for(R_xlen_t i(0); i < nsim; ++i) {
    samples[i] = ppjsdm::make_R_configuration(cpp_samples[i], types);
    Rcpp::as<Rcpp::List>(samples[i]).attr("number_points") = number_points[i];
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
                SEXP starting_configuration,
                R_xlen_t seed,
                R_xlen_t nthreads,
                bool debug) {
  using Configuration_type = std::vector<ppjsdm::Marked_point>;

  const auto cpp_window(ppjsdm::Window(window, mark_range));
  const ppjsdm::Truncated_exponential_family_model_over_window<std::vector<double>> exponential_model(cpp_window,
                                                                                                      Rcpp::as<std::vector<double>>(beta0),
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
    // TODO: Do not discard debug
    return rgibbs_helper<Configuration_type>(nthreads, seed, nsim, types, drop, exponential_model);
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
      return rgibbs_helper<Configuration_type>(nthreads,
                                               seed,
                                               nsim,
                                               types,
                                               drop,
                                               exponential_model,
                                               steps,
                                               debug,
                                               vector_starting_configuration,
                                               Rcpp::as<std::vector<R_xlen_t>>(only_simulate_these_types),
                                               vector_conditional_configuration);
    } else {
      return rgibbs_helper<Configuration_type>(nthreads,
                                               seed,
                                               nsim,
                                               types,
                                               drop,
                                               exponential_model,
                                               steps,
                                               debug,
                                               Rcpp::as<std::vector<R_xlen_t>>(only_simulate_these_types),
                                               vector_conditional_configuration);
    }
  }
}
