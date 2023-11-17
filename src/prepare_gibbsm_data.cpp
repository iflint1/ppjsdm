#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_manipulation.hpp"
#include "configuration/configuration_wrapper.hpp"
#include "configuration/get_number_points.hpp"
#include "configuration/make_R_configuration.hpp"

#include "point/point_manipulation.hpp"

#include "saturated_model/compute_dispersion_fitting.hpp"
#include "saturated_model/regression_vector.hpp"

#include "simulation/rbinomialpp_single.hpp"
#include "simulation/rstrat_single.hpp"

#include "utility/construct_if_missing.hpp"
#include "utility/gibbsm_helper_functions.hpp"
#include "utility/im_wrapper.hpp"
#include "utility/is_symmetric_matrix.hpp"
#include "utility/lightweight_matrix.hpp"
#include "utility/subset_to_nonNA.hpp"
#include "utility/timer.hpp"
#include "utility/window.hpp"

#include <algorithm> // std::remove_if, std::max_element, std::fill
#include <cmath> // std::log
#include <iterator> // std::next
#include <string> // std::string, std::to_string
#include <type_traits> // std::remove_cv_t
#include <vector> // std::vector

namespace detail {

// TODO: This function is copy/pasted from compute_vcov, and a similar version for numeric matrices
// was also used elsewhere; synchronise.
inline auto convert_list_of_boolean_matrices_to_vector(Rcpp::List boolean_matrices) {
  using conversion_t = ppjsdm::Lightweight_matrix<bool>;
  std::vector<conversion_t> converted_matrices(boolean_matrices.size());
  for(decltype(converted_matrices.size()) i(0); i < converted_matrices.size(); ++i) {
    converted_matrices[i] = conversion_t(Rcpp::as<Rcpp::LogicalMatrix>(boolean_matrices[i]));
  }
  return converted_matrices;
}

} // detail

template<typename Configuration, typename Vector>
auto make_dummy_points(const ppjsdm::Window& window,
                       const Vector& rho_times_volume,
                       std::string dummy_distribution) {
  Configuration D;
  if(dummy_distribution == std::string("binomial")) {
    D = ppjsdm::rbinomialpp_single<Configuration>(window, rho_times_volume);
  } else if(dummy_distribution == std::string("stratified")) {
    D = ppjsdm::rstratpp_single<Configuration>(window, rho_times_volume);
  } else {
    Rcpp::stop("Unexpected dummy_distribution parameter.");
  }

  return D;
}

template<typename computation_t, typename Vector>
auto make_rho_times_volume(int number_types,
                           const Vector& max_points_by_type,
                           R_xlen_t max_dummy,
                           R_xlen_t min_dummy,
                           double dummy_factor) {
  // This choice of intensity of dummy points rho is the guideline from the Baddeley et al. paper, see p. 8 therein.
  if(max_dummy == 0) {
    max_dummy = 1000000;
  }
  if(min_dummy == 0) {
    min_dummy = 500;
  }
  if(dummy_factor == 0.) {
    dummy_factor = 4.;
  }

  // Compute rho.
  using big_int_t = long long int;
  std::vector<big_int_t> rho_times_volume(number_types);

  for(typename decltype(rho_times_volume)::size_type i(0); i < rho_times_volume.size(); ++i) {
    const auto factor_times_max_points_in_type = static_cast<computation_t>(dummy_factor) * static_cast<computation_t>(max_points_by_type[i]);
    if(static_cast<big_int_t>(factor_times_max_points_in_type) < static_cast<big_int_t>(max_dummy)) {
      if(static_cast<big_int_t>(factor_times_max_points_in_type) > static_cast<big_int_t>(min_dummy)) {
        rho_times_volume[i] = static_cast<big_int_t>(factor_times_max_points_in_type);
      } else {
        rho_times_volume[i] = static_cast<big_int_t>(min_dummy);
      }
    } else {
      rho_times_volume[i] = static_cast<big_int_t>(max_dummy);
    }
  }
  return rho_times_volume;
}

template<typename Configuration, typename Dummy, typename FloatType, typename Vector, typename VectorBoolMatrices>
Rcpp::List prepare_gibbsm_data_helper(const std::vector<Configuration>& configuration_list,
                                      const Dummy& D,
                                      const Vector& rho_times_volume,
                                      const ppjsdm::Window& window,
                                      const ppjsdm::Im_list_wrapper& covariates,
                                      const std::vector<ppjsdm::Saturated_model<FloatType>>& dispersion_model,
                                      const ppjsdm::Saturated_model<FloatType>& medium_dispersion_model,
                                      const VectorBoolMatrices& estimate_alpha,
                                      Rcpp::LogicalMatrix estimate_gamma,
                                      int number_types,
                                      int nthreads,
                                      bool debug,
                                      Rcpp::CharacterVector type_names) {
  using size_t = ppjsdm::size_t<Configuration>;

  // Set up timer for debugging purposes
  ppjsdm::PreciseTimer timer{};
  if(debug) {
    Rcpp::Rcout << "Starting computation of the regression matrix...\n";
  }

  // Check if we actually need to compute alpha/gamma regressors
  bool need_to_compute_alpha(false), need_to_compute_gamma(false);
  for(decltype(estimate_gamma.nrow()) i(0); i < estimate_gamma.nrow(); ++i) {
    for(decltype(estimate_gamma.ncol()) j(i); j < estimate_gamma.ncol(); ++j) {
      for(decltype(estimate_alpha.size()) k(0); k < estimate_alpha.size(); ++k) {
        if(estimate_alpha[k](i, j)) {
          need_to_compute_alpha = true;
          break;
        }
      }
      if(estimate_gamma(i, j)) {
        need_to_compute_gamma = true;
      }
    }
    if(need_to_compute_alpha && need_to_compute_gamma) {
      break;
    }
  }

  // Remove points in D which generate NA values for one of the regressors.
  Dummy D_no_NA_covariates(ppjsdm::subset_to_nonNA(D, covariates));

  // Make shift vector
  Rcpp::NumericVector shift(Rcpp::no_init(number_types));
  const auto volume = static_cast<double>(window.volume());
  for(int j(0); j < number_types; ++j) {
    shift[j] = static_cast<double>(-std::log(static_cast<double>(rho_times_volume[j]) / volume));
  }

  // Precompute short and medium range dispersions
  using vec_dispersion_t = std::vector<decltype(ppjsdm::compute_dispersion_for_fitting(dispersion_model[0], 1, 1, configuration_list[0], D_no_NA_covariates))>;
  std::vector<vec_dispersion_t> dispersion_short{};
  for(decltype(dispersion_model.size()) j(0); j < dispersion_model.size(); ++j) {
    dispersion_short.emplace_back(vec_dispersion_t(configuration_list.size()));
  }
  vec_dispersion_t dispersion_medium(configuration_list.size());
  if(debug) {
    timer.set_current();
    Rcpp::Rcout << "Starting pre-computation of the dispersions to put into the regression matrix...\n";
  }
  for(decltype(configuration_list.size()) i(0); i < configuration_list.size(); ++i) {
    if(need_to_compute_alpha) {
      for(decltype(dispersion_short.size()) j(0); j < dispersion_short.size(); ++j) {
        dispersion_short[j][i] = ppjsdm::compute_dispersion_for_fitting(dispersion_model[j],
                                                                        number_types,
                                                                        nthreads,
                                                                        configuration_list[i], D_no_NA_covariates);
      }
    }
    if(need_to_compute_gamma) {
      dispersion_medium[i] = ppjsdm::compute_dispersion_for_fitting(medium_dispersion_model,
                                                                    number_types,
                                                                    nthreads,
                                                                    configuration_list[i], D_no_NA_covariates);
    }
  }
  if(debug) {
    Rcpp::Rcout << "Finished computing the dispersions. Time elapsed: " << timer.print_elapsed_time();
  }

  // Precompute how many of the points in the configuration we'll have to drop due to NA values on the covariates.
  size_t total_configuration_length(0);
  for(size_t i(0); i < configuration_list.size(); ++i) {
    total_configuration_length += ppjsdm::size(configuration_list[i]);
  }

  // Compute a few parameters necessary for the construction of regressors
  const auto total_points(total_configuration_length + D_no_NA_covariates.size() * configuration_list.size());
  const auto number_parameters_struct(ppjsdm::get_number_parameters(number_types,
                                                                    covariates.size(),
                                                                    estimate_alpha,
                                                                    estimate_gamma));
  const auto index_start_gamma(number_parameters_struct.index_start_gamma);
  const auto index_start_covariates(number_parameters_struct.index_start_covariates);
  const auto total_parameters(number_parameters_struct.total_parameters);

  // Default-initialise everything
  Rcpp::IntegerVector response(total_points);
  Rcpp::NumericVector x(total_points);
  Rcpp::NumericVector y(total_points);
  Rcpp::NumericVector mark(total_points);
  Rcpp::IntegerVector type(total_points);
  Rcpp::NumericMatrix regressors(total_points, total_parameters);

  // Fill the regressors, response, offset and shift with what we precomputed.
  size_t filling(0);
  for(size_t configuration_index(0); configuration_index < configuration_list.size(); ++configuration_index) {
    for(size_t point_index(0); point_index < ppjsdm::size(configuration_list[configuration_index]) + ppjsdm::size(D_no_NA_covariates); ++point_index) {
      const auto current_point(point_index < ppjsdm::size(configuration_list[configuration_index]) ?
                                 configuration_list[configuration_index][point_index] :
                                 D_no_NA_covariates[point_index - ppjsdm::size(configuration_list[configuration_index])]);
      // Check if the points generates any NA values, move on if so
      std::vector<decltype(covariates[0](current_point))> covariates_values(covariates.size());
      for(size_t covariate_index(0); covariate_index < static_cast<size_t>(covariates.size()); ++covariate_index) {
        const auto covariate(covariates[covariate_index](current_point));
        if(R_IsNA(covariate)) {
          Rcpp::warning("Some of the covariates had an NA value on points in the configuration. This should not happen, as the configurations have been subset to non-NA points.");
        }
        covariates_values[covariate_index] = covariate;
      }

      // fill response
      if(point_index < ppjsdm::size(configuration_list[configuration_index])) {
        response[filling] = 1;
      }

      // fill x, y, mark
      x[filling] = ppjsdm::get_x(current_point);
      y[filling] = ppjsdm::get_y(current_point);
      mark[filling] = ppjsdm::get_mark(current_point);

      // fill type
      const int type_index(ppjsdm::get_type(current_point));
      type[filling] = type_index + 1;

      // Fill regressors
      using current_dispersion_t = std::remove_cv_t<std::remove_reference_t<decltype(dispersion_short[0][0][0])>>;
      std::vector<current_dispersion_t> current_short{};
      current_dispersion_t current_medium;
      if(need_to_compute_alpha) {
        for(decltype(estimate_alpha.size()) k(0); k < estimate_alpha.size(); ++k) {
          current_short.emplace_back(dispersion_short[k][configuration_index][point_index]);
        }
      }
      if(need_to_compute_gamma) {
        current_medium = dispersion_medium[configuration_index][point_index];
      }
      const auto regression_vector(ppjsdm::make_regression_vector<double>(type_index,
                                                                          number_types,
                                                                          regressors.ncol(),
                                                                          index_start_covariates,
                                                                          index_start_gamma,
                                                                          covariates_values,
                                                                          current_short,
                                                                          current_medium,
                                                                          estimate_alpha,
                                                                          estimate_gamma));
      for(decltype(regression_vector.size()) k(0); k < regression_vector.size(); ++k) {
        regressors(filling, k) = regression_vector[k];
      }

      ++filling;
    }
  }

  const auto col_names(make_model_coloumn_names(covariates, estimate_alpha, estimate_gamma));
  Rcpp::colnames(regressors) = col_names;

  if(debug) {
    Rcpp::Rcout << "Finished computing the regression matrix. Time elapsed: " << timer.print_total_time();
  }

  return Rcpp::List::create(Rcpp::Named("response") = response,
                            Rcpp::Named("x") = x,
                            Rcpp::Named("y") = y,
                            Rcpp::Named("mark") = mark,
                            Rcpp::Named("type") = type,
                            Rcpp::Named("regressors") = regressors,
                            Rcpp::Named("shift") = shift,
                            Rcpp::Named("dummy") = ppjsdm::make_R_configuration(D_no_NA_covariates, type_names)
  );
}

// [[Rcpp::export]]
Rcpp::List prepare_gibbsm_data(Rcpp::List configuration_list,
                               SEXP window,
                               Rcpp::List covariates,
                               Rcpp::List model,
                               Rcpp::CharacterVector medium_range_model,
                               Rcpp::List short_range,
                               SEXP medium_range,
                               SEXP long_range,
                               R_xlen_t saturation,
                               Rcpp::NumericVector mark_range,
                               R_xlen_t max_dummy,
                               R_xlen_t min_dummy,
                               double dummy_factor,
                               Rcpp::List estimate_alpha,
                               Rcpp::LogicalMatrix estimate_gamma,
                               int nthreads,
                               bool debug,
                               std::string dummy_distribution,
                               Rcpp::CharacterVector type_names) {
  // TODO: A lot of code duplication between this function and the one below
  using computation_t = double;

  // Construct std::vector of configurations.
  std::vector<std::vector<ppjsdm::Marked_point>> vector_configurations(configuration_list.size());
  for(R_xlen_t i(0); i < configuration_list.size(); ++i) {
    const ppjsdm::Configuration_wrapper wrapped_configuration(Rcpp::as<Rcpp::List>(configuration_list[i]));
    const auto length_configuration(ppjsdm::size(wrapped_configuration));

    // Convert configurations to std::vector in order for parallelised version to work.
    vector_configurations[i] = std::vector<ppjsdm::Marked_point>(length_configuration);
    for(decltype(ppjsdm::size(wrapped_configuration)) j(0); j < length_configuration; ++j) {
      vector_configurations[i][j] = wrapped_configuration[j];
    }
  }

  const auto cpp_window(ppjsdm::Window(window, mark_range));

  auto max_points_by_type(ppjsdm::get_number_points(vector_configurations[0]));
  const auto number_types(max_points_by_type.size());

  using size_t = typename decltype(vector_configurations)::size_type;
  for(size_t i(0); i < vector_configurations.size(); ++i) {
    const auto current_points_by_type(ppjsdm::get_number_points(vector_configurations[i]));
    for(decltype(max_points_by_type.size()) j(0); j < number_types; ++j) {
      if(max_points_by_type[j] < current_points_by_type[j]) {
        max_points_by_type[j] = current_points_by_type[j];
      }
    }
  }
  using dispersion_t = decltype(ppjsdm::Saturated_model<double>(model[0], Rcpp::wrap(short_range[0]), saturation));
  std::vector<dispersion_t> dispersion{};
  for(decltype(short_range.size()) i(0); i < short_range.size(); ++i) {
    dispersion.emplace_back(ppjsdm::Saturated_model<double>(model[i], short_range[i], saturation));
  }
  const auto medium_range_dispersion(ppjsdm::Saturated_model<double>(medium_range_model, medium_range, long_range, saturation));

  // Compute rho.
  const auto rho_times_volume(make_rho_times_volume<computation_t>(number_types,
                                                                   max_points_by_type,
                                                                   max_dummy,
                                                                   min_dummy,
                                                                   dummy_factor));

  // Sample the dummy points D.
  auto D(make_dummy_points<std::vector<ppjsdm::Marked_point>>(cpp_window,
                                                              rho_times_volume,
                                                              dummy_distribution));

  return prepare_gibbsm_data_helper(vector_configurations,
                                    D,
                                    rho_times_volume,
                                    cpp_window,
                                    ppjsdm::Im_list_wrapper(covariates),
                                    dispersion,
                                    medium_range_dispersion,
                                    detail::convert_list_of_boolean_matrices_to_vector(estimate_alpha),
                                    estimate_gamma,
                                    number_types,
                                    nthreads,
                                    debug,
                                    type_names);
}

// [[Rcpp::export]]
Rcpp::List prepare_gibbsm_data_with_dummy(Rcpp::List configuration_list,
                                          SEXP dummy,
                                          SEXP window,
                                          Rcpp::List covariates,
                                          Rcpp::List model,
                                          Rcpp::CharacterVector medium_range_model,
                                          Rcpp::List short_range,
                                          SEXP medium_range,
                                          SEXP long_range,
                                          R_xlen_t saturation,
                                          Rcpp::NumericVector mark_range,
                                          Rcpp::List estimate_alpha,
                                          Rcpp::LogicalMatrix estimate_gamma,
                                          int nthreads,
                                          bool debug,
                                          Rcpp::CharacterVector type_names) {
  // Construct std::vector of configurations.
  std::vector<std::vector<ppjsdm::Marked_point>> vector_configurations(configuration_list.size());
  for(R_xlen_t i(0); i < configuration_list.size(); ++i) {
    const ppjsdm::Configuration_wrapper wrapped_configuration(Rcpp::as<Rcpp::List>(configuration_list[i]));
    const auto length_configuration(ppjsdm::size(wrapped_configuration));

    // Convert configurations to std::vector in order for parallelised version to work.
    vector_configurations[i] = std::vector<ppjsdm::Marked_point>(length_configuration);
    for(decltype(ppjsdm::size(wrapped_configuration)) j(0); j < length_configuration; ++j) {
      vector_configurations[i][j] = wrapped_configuration[j];
    }
  }

  // Construct std::vector of dummy
  const ppjsdm::Configuration_wrapper wrapped_dummy(Rcpp::as<Rcpp::List>(dummy));
  const auto length_dummy(ppjsdm::size(wrapped_dummy));

  std::vector<ppjsdm::Marked_point> vector_dummy(length_dummy);
  for(decltype(ppjsdm::size(wrapped_dummy)) j(0); j < length_dummy; ++j) {
    vector_dummy[j] = wrapped_dummy[j];
  }

  const auto cpp_window(ppjsdm::Window(window, mark_range));

  auto max_points_by_type(ppjsdm::get_number_points(vector_configurations[0]));
  const auto number_types(max_points_by_type.size());

  using size_t = typename decltype(vector_configurations)::size_type;
  for(size_t i(0); i < vector_configurations.size(); ++i) {
    const auto current_points_by_type(ppjsdm::get_number_points(vector_configurations[i]));
    for(decltype(max_points_by_type.size()) j(0); j < number_types; ++j) {
      if(max_points_by_type[j] < current_points_by_type[j]) {
        max_points_by_type[j] = current_points_by_type[j];
      }
    }
  }
  using dispersion_t = decltype(ppjsdm::Saturated_model<double>(model[0], Rcpp::wrap(short_range[0]), saturation));
  std::vector<dispersion_t> dispersion{};
  for(decltype(short_range.size()) i(0); i < short_range.size(); ++i) {
    dispersion.emplace_back(ppjsdm::Saturated_model<double>(model[i], short_range[i], saturation));
  }
  const auto medium_range_dispersion(ppjsdm::Saturated_model<double>(medium_range_model, medium_range, long_range, saturation));

  return prepare_gibbsm_data_helper(vector_configurations,
                                    vector_dummy,
                                    ppjsdm::get_number_points(vector_dummy),
                                    cpp_window,
                                    ppjsdm::Im_list_wrapper(covariates),
                                    dispersion,
                                    medium_range_dispersion,
                                    detail::convert_list_of_boolean_matrices_to_vector(estimate_alpha),
                                    estimate_gamma,
                                    number_types,
                                    nthreads,
                                    debug,
                                    type_names);
}
