#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_manipulation.hpp"
#include "configuration/configuration_wrapper.hpp"
#include "configuration/get_number_points.hpp"

#include "point/point_manipulation.hpp"

#include "saturated_model/saturated_model.hpp"

#include "simulation/rppp_single.hpp"

#include "utility/construct_if_missing.hpp"
#include "utility/gibbsm_helper_functions.hpp"
#include "utility/im_wrapper.hpp"
#include "utility/is_symmetric_matrix.hpp"
#include "utility/lightweight_matrix.hpp"
#include "utility/window.hpp"

#include <algorithm> // std::remove_if, std::max_element, std::fill
#include <cmath> // std::log
#include <iterator> // std::next
#include <string> // std::string, std::to_string
#include <sstream> // std::stringstream
#include <type_traits> // std::remove_cv_t
#include <vector> // std::vector

#ifdef _OPENMP
#include <omp.h>
#endif

inline void add_to_formula(std::string& formula, Rcpp::CharacterVector names) {
  for(const auto& name: names) {
    formula += std::string(" + ") + Rcpp::as<std::string>(name);
  }
}

template<typename Configuration, typename Vector>
Rcpp::List prepare_gibbsm_data_helper(const std::vector<Configuration>& configuration_list,
                                      const ppjsdm::Window& window,
                                      const ppjsdm::Im_list_wrapper& covariates,
                                      const ppjsdm::Saturated_model& dispersion_model,
                                      const ppjsdm::Saturated_model& medium_dispersion_model,
                                      const Vector& max_points_by_type,
                                      R_xlen_t max_dummy,
                                      double dummy_factor,
                                      Rcpp::LogicalMatrix estimate_alpha,
                                      Rcpp::LogicalMatrix estimate_gamma,
                                      int number_types,
                                      int nthreads) {
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif

  // Check if we actually need to compute alpha/gamma regressors
  bool need_to_compute_alpha(false);
  bool need_to_compute_gamma(false);
  for(R_xlen_t i(0); i < estimate_alpha.nrow(); ++i) {
    for(R_xlen_t j(0); j < estimate_alpha.ncol(); ++j) {
      if(estimate_alpha(i, j)) {
        need_to_compute_alpha = true;
      }
      if(estimate_gamma(i, j)) {
        need_to_compute_gamma = true;
      }
    }
  }

  using size_t = ppjsdm::size_t<Configuration>;

  // This choice of intensity of dummy points rho is the guideline from the Baddeley et al. paper, see p. 8 therein.
  if(max_dummy == 0) {
    max_dummy = 500;
  }
  if(dummy_factor == 0.) {
    dummy_factor = 4.;
  }

  // Sample the dummy points D.
  std::vector<double> rho_times_volume(number_types);
  size_t length_D(0);
  for(typename decltype(rho_times_volume)::size_type i(0); i < rho_times_volume.size(); ++i) {
    const auto factor_times_max_points_in_any_type(dummy_factor * static_cast<double>(max_points_by_type[i]));
    if(factor_times_max_points_in_any_type < static_cast<double>(max_dummy)) {
      rho_times_volume[i] = factor_times_max_points_in_any_type;
      length_D += static_cast<size_t>(factor_times_max_points_in_any_type);
    } else {
      rho_times_volume[i] = max_dummy;
      length_D += max_dummy;
    }
  }
  auto D(ppjsdm::rbinomialpp_single<std::vector<ppjsdm::Marked_point>>(window, rho_times_volume, number_types, length_D));

  // Remove points in D which generate NA values for one of the regressors.
  D.erase(std::remove_if(D.begin(), D.end(),
      [&covariates](const auto& point) {
        for(std::remove_cv_t<decltype(covariates.size())> covariate_index(0); covariate_index < covariates.size(); ++covariate_index) {
          const auto covariate(covariates[covariate_index](point));
          if(R_IsNA(covariate)) {
            return true;
          }
        }
        return false;
      }), D.end());
  length_D = D.size();

  // Make shift vector
  Rcpp::NumericVector shift(Rcpp::no_init(number_types));
  const auto volume(window.volume());
  for(int j(0); j < number_types; ++j) {
    shift[j] = -std::log(static_cast<double>(rho_times_volume[j]) / volume);
  }

  // Precompute short and medium range dispersions
  std::vector<std::vector<std::vector<double>>> dispersion_short(configuration_list.size());
  std::vector<std::vector<std::vector<double>>> dispersion_medium(configuration_list.size());
  for(size_t configuration_index(0); configuration_index < configuration_list.size(); ++configuration_index) {
    if(need_to_compute_alpha) {
      dispersion_short[configuration_index] = ppjsdm::compute_dispersion_index(dispersion_model,
                                                                               number_types,
                                                                               configuration_list[configuration_index], D);
    }
    if(need_to_compute_gamma) {
      dispersion_medium[configuration_index] = ppjsdm::compute_dispersion_index(medium_dispersion_model,
                                                                                number_types,
                                                                                configuration_list[configuration_index], D);
    }
  }

  // Precompute how many of the points in the configuration we'll have to drop due to NA values on the covariates.
  unsigned long long int total_configuration_length(0);
  size_t total_removed_points(0);
  for(size_t configuration_index(0); configuration_index < configuration_list.size(); ++configuration_index) {
    total_configuration_length += ppjsdm::size(configuration_list[configuration_index]);
    for(size_t point_index(0); point_index < ppjsdm::size(configuration_list[configuration_index]); ++point_index) {
      for(std::remove_cv_t<decltype(covariates.size())> covariate_index(0); covariate_index < covariates.size(); ++covariate_index) {
        const auto covariate(covariates[covariate_index](configuration_list[configuration_index][point_index]));
        if(R_IsNA(covariate)) {
          ++total_removed_points;
          break;
        }
      }
    }
  }

  // Send a warning if we had to drop any
  if(total_removed_points > 0) {
    std::stringstream ss;
    ss << "Some of the covariates had an NA value on " << total_removed_points << " of the configuration points; these have been removed.";
    Rcpp::warning(ss.str());
  }

  // Compute a few parameters necessary for the construction of regressors
  const auto total_points(total_configuration_length - total_removed_points + length_D * configuration_list.size());
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
  Rcpp::NumericVector rho_offset(total_points);
  Rcpp::NumericMatrix regressors(total_points, total_parameters);

  // Fill the regressors, response, offset and shift with what we precomputed.
  size_t filling(0);
  for(size_t configuration_index(0); configuration_index < configuration_list.size(); ++configuration_index) {
    for(size_t point_index(0); point_index < ppjsdm::size(configuration_list[configuration_index]) + ppjsdm::size(D); ++point_index) {
      const auto current_point(point_index < ppjsdm::size(configuration_list[configuration_index]) ?
                                 configuration_list[configuration_index][point_index] :
                                 D[point_index - ppjsdm::size(configuration_list[configuration_index])]);
      // Check if the points generates any NA values, move on if so
      std::vector<double> covariates_values(covariates.size());
      bool found_na_point(false);
      for(size_t covariate_index(0); covariate_index < static_cast<size_t>(covariates.size()); ++covariate_index) {
        const auto covariate(covariates[covariate_index](current_point));
        if(R_IsNA(covariate)) {
          found_na_point = true;
          break;
        }
        covariates_values[covariate_index] = covariate;
      }
      if(found_na_point) {
        continue;
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

      // fill rho
      rho_offset[filling] = -std::log(static_cast<double>(rho_times_volume[type_index]) / volume);

      // fill log_lambda
      regressors(filling, type_index) = 1.;

      // fill covariates
      for(size_t covariate_index(0); covariate_index < static_cast<size_t>(covariates.size()); ++covariate_index) {
        regressors(filling, index_start_covariates + covariate_index * number_types + type_index) = covariates_values[covariate_index];
      }

      // fill alpha & gamma
      size_t index_alpha(0);
      size_t index_gamma(0);
      for(int j(0); j < number_types; ++j) {
        for(int k(j); k < number_types; ++k) {
          if(estimate_alpha(j, k)) {
            if(j == type_index) {
              regressors(filling, number_types + index_alpha) = dispersion_short[configuration_index][point_index][k];
            } else if(k == type_index) {
              regressors(filling, number_types + index_alpha) = dispersion_short[configuration_index][point_index][j];
            }
            ++index_alpha;
          }
          if(estimate_gamma(j, k)) {
            if(j == type_index) {
              regressors(filling, index_start_gamma + index_gamma) = dispersion_medium[configuration_index][point_index][k];
            } else if(k == type_index) {
              regressors(filling, index_start_gamma + index_gamma) = dispersion_medium[configuration_index][point_index][j];
            }
            ++index_gamma;
          }
        }
      }
      ++filling;
    }
  }

  const auto col_names(make_model_coloumn_names(covariates, number_types, estimate_alpha, estimate_gamma));
  Rcpp::colnames(regressors) = col_names;

  return Rcpp::List::create(Rcpp::Named("response") = response,
                            Rcpp::Named("x") = x,
                            Rcpp::Named("y") = y,
                            Rcpp::Named("mark") = mark,
                            Rcpp::Named("type") = type,
                            Rcpp::Named("offset") = rho_offset,
                            Rcpp::Named("regressors") = regressors,
                            Rcpp::Named("shift") = shift
  );
}

// [[Rcpp::export]]
Rcpp::List prepare_gibbsm_data(Rcpp::List configuration_list,
                               SEXP window,
                               Rcpp::List covariates,
                               Rcpp::CharacterVector model,
                               Rcpp::CharacterVector medium_range_model,
                               SEXP short_range,
                               SEXP medium_range,
                               SEXP long_range,
                               R_xlen_t saturation,
                               Rcpp::NumericVector mark_range,
                               R_xlen_t max_dummy,
                               double dummy_factor,
                               Rcpp::LogicalMatrix estimate_alpha,
                               Rcpp::LogicalMatrix estimate_gamma,
                               int nthreads) {
  // Construct std::vector of configurations.
  std::vector<std::vector<ppjsdm::Marked_point>> vector_configurations(configuration_list.size());
  for(R_xlen_t i(0); i < configuration_list.size(); ++i) {
    const ppjsdm::Configuration_wrapper wrapped_configuration(Rcpp::wrap(configuration_list[i]));
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
    const auto new_points_by_type(ppjsdm::get_number_points(vector_configurations[i]));
    for(size_t j(0); j < max_points_by_type.size(); ++j) {
      if(max_points_by_type[j] < new_points_by_type[j]) {
        max_points_by_type[j] = new_points_by_type[j];
      }
    }
  }
  const auto dispersion(ppjsdm::Saturated_model(model, short_range, saturation));
  const auto medium_range_dispersion(ppjsdm::Saturated_model(medium_range_model, medium_range, long_range, saturation));
  return prepare_gibbsm_data_helper(vector_configurations, cpp_window, ppjsdm::Im_list_wrapper(covariates), dispersion, medium_range_dispersion, max_points_by_type, max_dummy, dummy_factor, estimate_alpha, estimate_gamma, number_types, nthreads);
}
