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

  // Sample the dummy points D.
  // This choice of rho is the guideline from the Baddeley et al. paper, see p. 8 therein.
  if(max_dummy == 0) {
    max_dummy = 500;
  }
  if(dummy_factor == 0.) {
    dummy_factor = 4.;
  }
  std::vector<double> rho_times_volume(number_types);
  size_t length_D(0);
  for(typename decltype(rho_times_volume)::size_type i(0); i < rho_times_volume.size(); ++i) {
    const int factor_times_max_points_in_any_type(dummy_factor * max_points_by_type[i]);
    if(factor_times_max_points_in_any_type < max_dummy) {
      rho_times_volume[i] = factor_times_max_points_in_any_type;
      length_D += factor_times_max_points_in_any_type;
    } else {
      rho_times_volume[i] = max_dummy;
      length_D += max_dummy;
    }
  }

  const auto D(ppjsdm::rbinomialpp_single<std::vector<ppjsdm::Marked_point>>(window, rho_times_volume, number_types, length_D));

  // Make shift vector
  Rcpp::NumericVector shift(Rcpp::no_init(number_types));
  const auto volume(window.volume());
  for(int j(0); j < number_types; ++j) {
    shift[j] = -std::log(static_cast<double>(rho_times_volume[j]) / volume);
  }

  // Precompute dispersion and other computation-intensive things.
  struct Result_type {
    Result_type(): is_in_configuration(), type(), dispersion(), medium_dispersion(), covariates() {}

    Result_type(bool i, int t, std::vector<double> v, std::vector<double> w, std::vector<double> x, double new_x, double new_y, double m):
      is_in_configuration(i), type(t), dispersion(std::move(v)), medium_dispersion(std::move(w)), covariates(std::move(x)), x(new_x), y(new_y), mark(m) {}

    bool is_in_configuration;
    int type;
    std::vector<double> dispersion;
    std::vector<double> medium_dispersion;
    std::vector<double> covariates;
    double x;
    double y;
    double mark;
  };
  std::vector<Result_type> precomputed_results;
  unsigned long long int total_configuration_length(0);
  for(size_t i(0); i < configuration_list.size(); ++i) {
    total_configuration_length += ppjsdm::size(configuration_list[i]);
  }
  precomputed_results.reserve(total_configuration_length + length_D * configuration_list.size());

  const size_t covariates_length(covariates.size());
#pragma omp parallel
  {
  std::vector<Result_type> results_private;
  results_private.reserve(total_configuration_length + length_D * configuration_list.size());
#pragma omp for nowait
  for(size_t i = 0; i < total_configuration_length; ++i) {
    size_t configuration_index(0);
    unsigned long long int previous_count(0);
    while(i >= previous_count + ppjsdm::size(configuration_list[configuration_index])) {
      previous_count += ppjsdm::size(configuration_list[configuration_index]);
      ++configuration_index;
    }
    // configuration_index contains the index of the configuration we're at in the loop.
    const auto point_index(i - previous_count);
    std::vector<double> d;
    if(need_to_compute_alpha) {
      d = ppjsdm::compute_dispersion(dispersion_model, configuration_list[configuration_index][point_index], number_types, configuration_list[configuration_index]);
    }
    std::vector<double> e;
    if(need_to_compute_gamma) {
      e = ppjsdm::compute_dispersion(medium_dispersion_model, configuration_list[configuration_index][point_index], number_types, configuration_list[configuration_index]);
    }
    std::vector<double> cov(covariates_length);
    for(size_t k(0); k < covariates_length; ++k) {
      const auto covariate(covariates[k](configuration_list[configuration_index][point_index]));
      if(R_IsNA(covariate)) {
        Rcpp::stop("One of the covariates' value is NA on one of the locations in the dataset.");
      }
      cov[k] = covariate;
    }
    results_private.emplace_back(true, ppjsdm::get_type(configuration_list[configuration_index][point_index]), std::move(d), std::move(e), std::move(cov), ppjsdm::get_x(configuration_list[configuration_index][point_index]), ppjsdm::get_y(configuration_list[configuration_index][point_index]), ppjsdm::get_mark(configuration_list[configuration_index][point_index]));
  }
#pragma omp for nowait
  for(size_t i = 0; i < length_D; ++i) {
    bool one_of_covariates_na(false);
    std::vector<double> cov(covariates_length);
    for(size_t k(0); k < covariates_length; ++k) {
      const auto covariate(covariates[k](D[i]));
      if(R_IsNA(covariate)) {
        one_of_covariates_na = true;
        break;
      }
      cov[k] = covariate;
    }
    // Get rid of locations with an NA value for one of the covariates.
    if(one_of_covariates_na) {
      continue;
    }
    for(size_t j(0); j < configuration_list.size(); ++j) {
      std::vector<double> d;
      if(need_to_compute_alpha) {
        d = ppjsdm::compute_dispersion(dispersion_model, D[i], number_types, configuration_list[j]);
      }
      std::vector<double> e;
      if(need_to_compute_gamma) {
        e = ppjsdm::compute_dispersion(medium_dispersion_model, D[i], number_types, configuration_list[j]);
      }
      results_private.emplace_back(false, ppjsdm::get_type(D[i]), std::move(d), std::move(e), cov, ppjsdm::get_x(D[i]), ppjsdm::get_y(D[i]), ppjsdm::get_mark(D[i]));
    }
  }
#pragma omp critical
  precomputed_results.insert(precomputed_results.end(), results_private.begin(), results_private.end());
  }

  //const size_t number_traits(traits.size());
  const auto total_points(precomputed_results.size());
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
  for(size_t i(0); i < total_points; ++i) {
    if(precomputed_results[i].is_in_configuration) {
      response[i] = 1;
    }
    x[i] = precomputed_results[i].x;
    y[i] = precomputed_results[i].y;
    mark[i] = precomputed_results[i].mark;

    const int type_index(precomputed_results[i].type);
    type[i] = type_index + 1;

    rho_offset[i] = -std::log(static_cast<double>(rho_times_volume[type_index]) / volume);

    // fill traits
    // short_range_traits_input(i, 0) = precomputed_results[i].dispersion[type_index];
    // medium_range_traits_input(i, 0) = precomputed_results[i].medium_dispersion[type_index];
    // for(int j(0); j < number_types; ++j) {
    //   if(j != type_index) {
    //     short_range_joint_traits_input(i, 0) += precomputed_results[i].dispersion[j];
    //     medium_range_joint_traits_input(i, 0) += precomputed_results[i].medium_dispersion[j];
    //   }
    // }
    //
    // for(size_t k(0); k < number_traits; ++k) {
    //   short_range_traits_input(i, k + 1) = Rcpp::as<Rcpp::NumericVector>(traits[k])[type_index] * precomputed_results[i].dispersion[type_index];
    //   medium_range_traits_input(i, k + 1) = Rcpp::as<Rcpp::NumericVector>(traits[k])[type_index] * precomputed_results[i].medium_dispersion[type_index];
    //   for(int j(0); j < number_types; ++j) {
    //     if(j != type_index) {
    //       const auto delta(std::abs(Rcpp::as<Rcpp::NumericVector>(traits[k])[j] - Rcpp::as<Rcpp::NumericVector>(traits[k])[type_index]));
    //       short_range_joint_traits_input(i, k + 1) += delta * precomputed_results[i].dispersion[j];
    //       medium_range_joint_traits_input(i, k + 1) += delta * precomputed_results[i].medium_dispersion[j];
    //     }
    //   }
    // }

    size_t index_alpha(0);
    size_t index_gamma(0);
    for(int j(0); j < number_types; ++j) {
      if(j == type_index) {
        // fill log_lambda
        regressors(i, j) = 1.;

        // fill covariates
        for(size_t k(0); k < covariates_length; ++k) {
          regressors(i, index_start_covariates + k * number_types + j) = precomputed_results[i].covariates[k];
        }

        // fill alpha & gamma
        for(int k(j); k < number_types; ++k) {
          if(estimate_alpha(j, k)) {
            regressors(i, number_types + index_alpha) = precomputed_results[i].dispersion[k];
            ++index_alpha;
          }
          if(estimate_gamma(j, k)) {
            regressors(i, index_start_gamma + index_gamma) = precomputed_results[i].medium_dispersion[k];
            ++index_gamma;
          }
        }
      } else {
        // fill alpha & gamma
        for(int k(j); k < number_types; ++k) {
          if(estimate_alpha(j, k)) {
            if(k == type_index) {
              regressors(i, number_types + index_alpha) = precomputed_results[i].dispersion[j];
            }
            ++index_alpha;
          }
          if(estimate_gamma(j, k)) {
            if(k == type_index) {
              regressors(i, index_start_gamma + index_gamma) = precomputed_results[i].medium_dispersion[j];
            }
            ++index_gamma;
          }
        }
      }
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
  //auto max_points_in_any_type(*std::max_element(max_points_by_type.begin(), max_points_by_type.end()));

  using size_t = typename decltype(vector_configurations)::size_type;
  for(size_t i(0); i < vector_configurations.size(); ++i) {
    const auto new_points_by_type(ppjsdm::get_number_points(vector_configurations[i]));
    for(size_t j(0); j < max_points_by_type.size(); ++j) {
      // if(max_points_in_any_type < new_points_by_type[j]) {
      //   max_points_in_any_type = new_points_by_type[j];
      // }
      if(max_points_by_type[j] < new_points_by_type[j]) {
        max_points_by_type[j] = new_points_by_type[j];
      }
    }
  }
  const auto dispersion(ppjsdm::Saturated_model(model, short_range, saturation));
  const auto medium_range_dispersion(ppjsdm::Saturated_model(medium_range_model, medium_range, long_range, saturation));
  return prepare_gibbsm_data_helper(vector_configurations, cpp_window, ppjsdm::Im_list_wrapper(covariates), dispersion, medium_range_dispersion, max_points_by_type, max_dummy, dummy_factor, estimate_alpha, estimate_gamma, number_types, nthreads);
}
