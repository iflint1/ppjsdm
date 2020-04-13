#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_manipulation.hpp"
#include "configuration/configuration_wrapper.hpp"
#include "configuration/get_number_points.hpp"

#include "point/point_manipulation.hpp"

#include "saturated_model/saturated_model.hpp"

#include "simulation/rbinomialpp_single.hpp"

#include "utility/construct_if_missing.hpp"
#include "utility/im_wrapper.hpp"
#include "utility/is_symmetric_matrix.hpp"
#include "utility/lightweight_matrix.hpp"
#include "utility/window_utilities.hpp"

#include <algorithm> // std::remove_if
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

template<bool Approximate, typename Configuration, typename DispersionModel, typename MediumDispersionModel, typename Vector>
Rcpp::List prepare_gibbsm_data_helper(const std::vector<Configuration>& configuration_list, const ppjsdm::Window& window, const ppjsdm::Im_list_wrapper& covariates, Rcpp::List traits, const DispersionModel& dispersion_model, const MediumDispersionModel& medium_dispersion_model, const Vector& max_points_by_type) {
  using size_t = ppjsdm::size_t<Configuration>;

  // Sample the dummy points D.
  // This choice of rho is the guideline from the Baddeley et al. paper, see p. 8 therein.
  Vector rho_times_volume(max_points_by_type);
  size_t length_D(0);
  for(auto& n: rho_times_volume) {
    const auto mult_by_four(n * 4);
    n = mult_by_four < 500 ? 500 : mult_by_four;
    length_D += n;
  }
  const size_t number_types(max_points_by_type.size());
  auto D(ppjsdm::rbinomialpp_single<std::vector<ppjsdm::Marked_point>>(window, rho_times_volume, number_types, length_D));

  // Make shift vector
  Rcpp::NumericVector shift(Rcpp::no_init(number_types));
  const auto volume(window.volume());
  for(size_t j(0); j < number_types; ++j) {
    shift[j] = -std::log(static_cast<double>(rho_times_volume[j]) / volume);
  }

  // Precompute dispersion and other computation-intensive things.
  struct Result_type {
    Result_type(): is_in_configuration(), type(), dispersion(), medium_dispersion(), covariates() {}

    Result_type(bool i, int t, std::vector<double> v, std::vector<double> w, std::vector<double> x):
      is_in_configuration(i), type(t), dispersion(std::move(v)), medium_dispersion(std::move(w)), covariates(std::move(x)) {}

    bool is_in_configuration;
    int type;
    std::vector<double> dispersion;
    std::vector<double> medium_dispersion;
    std::vector<double> covariates;
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
    // TODO: Avoidable?
    Configuration configuration_copy(configuration_list[configuration_index]);
    ppjsdm::remove_point_by_iterator(configuration_copy, std::next(configuration_copy.begin(), point_index));
    const auto d(ppjsdm::compute_dispersion<Approximate>(dispersion_model, configuration_list[configuration_index][point_index], number_types, configuration_copy));
    const auto e(ppjsdm::compute_dispersion<Approximate>(medium_dispersion_model, configuration_list[configuration_index][point_index], number_types, configuration_copy));
    std::vector<double> cov(covariates_length);
    for(size_t k(0); k < covariates_length; ++k) {
      const auto covariate(covariates[k](configuration_list[configuration_index][point_index]));
      if(R_IsNA(covariate)) {
        Rcpp::stop("One of the covariates' value is NA on one of the locations in the dataset.");
      }
      cov[k] = covariate;
    }
    results_private.emplace_back(true, ppjsdm::get_type(configuration_list[configuration_index][point_index]), std::move(d), std::move(e), std::move(cov));
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
      const auto d(ppjsdm::compute_dispersion<Approximate>(dispersion_model, D[i], number_types, configuration_list[j]));
      const auto e(ppjsdm::compute_dispersion<Approximate>(medium_dispersion_model, D[i], number_types, configuration_list[j]));
      results_private.emplace_back(false, ppjsdm::get_type(D[i]), std::move(d), std::move(e), cov);
    }
  }
#pragma omp critical
  precomputed_results.insert(precomputed_results.end(), results_private.begin(), results_private.end());
  }

  const auto total_points(precomputed_results.size());
  const size_t number_traits(traits.size());

  // Default-initialise the data
  std::vector<int> response(total_points);
  std::vector<double> rho_offset(total_points);
  ppjsdm::Lightweight_matrix<double> log_lambda(total_points, number_types);
  ppjsdm::Lightweight_matrix<double> alpha_input(total_points, number_types * (number_types + 1) / 2);
  ppjsdm::Lightweight_matrix<double> gamma_input(total_points, number_types * (number_types + 1) / 2);
  ppjsdm::Lightweight_matrix<double> covariates_input(total_points, covariates_length * number_types);
  ppjsdm::Lightweight_matrix<double> short_range_traits_input(total_points, 1 + number_traits);
  ppjsdm::Lightweight_matrix<double> short_range_joint_traits_input(total_points, 1 + number_traits);
  ppjsdm::Lightweight_matrix<double> medium_range_traits_input(total_points, 1 + number_traits);
  ppjsdm::Lightweight_matrix<double> medium_range_joint_traits_input(total_points, 1 + number_traits);


  // Fill the regressors, response, offset and shift with what we precomputed.
  for(size_t i(0); i < total_points; ++i) {
    response[i] = precomputed_results[i].is_in_configuration ? 1 : 0;
    const size_t type_index(precomputed_results[i].type);

    rho_offset[i] = -std::log(static_cast<double>(rho_times_volume[type_index]) / volume);

    // fill traits
    short_range_traits_input(i, 0) = precomputed_results[i].dispersion[type_index];
    medium_range_traits_input(i, 0) = precomputed_results[i].medium_dispersion[type_index];
    // TODO: Initialisation probably not needed, just testing some things out.
    short_range_joint_traits_input(i, 0) = 0;
    medium_range_joint_traits_input(i, 0) = 0;
    for(size_t j(0); j < number_types; ++j) {
      if(j != type_index) {
        short_range_joint_traits_input(i, 0) += precomputed_results[i].dispersion[j];
        medium_range_joint_traits_input(i, 0) += precomputed_results[i].medium_dispersion[j];
      }
    }

    for(size_t k(0); k < number_traits; ++k) {
      short_range_traits_input(i, k + 1) = Rcpp::as<Rcpp::NumericVector>(traits[k])[type_index] * precomputed_results[i].dispersion[type_index];
      medium_range_traits_input(i, k + 1) = Rcpp::as<Rcpp::NumericVector>(traits[k])[type_index] * precomputed_results[i].medium_dispersion[type_index];
      // TODO: Initialisation probably not needed, just testing some things out.
      short_range_joint_traits_input(i, k + 1) = 0;
      medium_range_joint_traits_input(i, k + 1) = 0;
      for(size_t j(0); j < number_types; ++j) {
        if(j != type_index) {
          const auto delta(std::abs(Rcpp::as<Rcpp::NumericVector>(traits[k])[j] - Rcpp::as<Rcpp::NumericVector>(traits[k])[type_index]));
          short_range_joint_traits_input(i, k + 1) += delta * precomputed_results[i].dispersion[j];
          medium_range_joint_traits_input(i, k + 1) += delta * precomputed_results[i].medium_dispersion[j];
        }
      }
    }

    // TODO: index or formula here and in other vectors?
    size_t index(0);
    for(size_t j(0); j < number_types; ++j) {
      if(j == type_index) {
        // fill log_lambda
        log_lambda(i, j) = 1;

        // fill covariates
        for(size_t k(0); k < covariates_length; ++k) {
          covariates_input(i, k * number_types + j) = precomputed_results[i].covariates[k];
        }

        // fill alpha
        for(size_t k(j); k < number_types; ++k) {
          alpha_input(i, index) = precomputed_results[i].dispersion[k];
          gamma_input(i, index++) = precomputed_results[i].medium_dispersion[k];
        }
      } else {
        // fill log_lambda
        log_lambda(i, j) = 0;

        // fill covariates
        for(size_t k(0); k < covariates_length; ++k) {
          covariates_input(i, k * number_types + j) = 0;
        }

        // fill alpha
        for(size_t k(j); k < number_types; ++k) {
          if(k == type_index) {
            alpha_input(i, index) = precomputed_results[i].dispersion[j];
            gamma_input(i, index++) = precomputed_results[i].medium_dispersion[j];
          } else {
            // TODO: Not necessary since zero-initialized.
            alpha_input(i, index) = 0;
            gamma_input(i, index++) = 0;
          }
        }
      }
    }
  }

  // Convert response and rho_offset to Rcpp objects.
  Rcpp::IntegerMatrix response_rcpp(Rcpp::no_init(response.size(), 1));
  Rcpp::NumericMatrix rho_offset_rcpp(Rcpp::no_init(response.size(), 1));
  for(R_xlen_t i(0); i < static_cast<R_xlen_t>(response.size()); ++i) {
    response_rcpp(i, 0) = response[i];
    rho_offset_rcpp(i, 0) = rho_offset[i];
  }

  // Set names.
  Rcpp::CharacterVector alpha_names(Rcpp::no_init(alpha_input.ncol()));
  Rcpp::CharacterVector gamma_names(Rcpp::no_init(alpha_input.ncol()));

  size_t index(0);
  for(size_t j(0); j < number_types; ++j) {
    for(size_t k(j); k < number_types; ++k) {
      alpha_names[index] = std::string("alpha_") + std::to_string(j + 1) + std::string("_") + std::to_string(k + 1);
      gamma_names[index++] = std::string("gamma_") + std::to_string(j + 1) + std::string("_") + std::to_string(k + 1);
    }
  }

  // TODO: Better names
  Rcpp::CharacterVector a_names(Rcpp::no_init(short_range_traits_input.ncol()));
  Rcpp::CharacterVector b_names(Rcpp::no_init(medium_range_traits_input.ncol()));
  Rcpp::CharacterVector c_names(Rcpp::no_init(short_range_joint_traits_input.ncol()));
  Rcpp::CharacterVector d_names(Rcpp::no_init(medium_range_joint_traits_input.ncol()));
  for(size_t i(0); i < short_range_traits_input.ncol(); ++i) {
    a_names[i] = std::string("short_range_direct_") + (i == 0 ? std::string("intercept") : std::to_string(i));
    b_names[i] = std::string("medium_range_direct_") + (i == 0 ? std::string("intercept") : std::to_string(i));
    c_names[i] = std::string("short_range_joint_") + (i == 0 ? std::string("intercept") : std::to_string(i));
    d_names[i] = std::string("medium_range_joint_") + (i == 0 ? std::string("intercept") : std::to_string(i));
  }

  Rcpp::CharacterVector covariates_input_names(Rcpp::no_init(covariates_length * number_types));
  if(covariates_length > 0) {
    const auto covariates_names(covariates.names());
    for(size_t i(0); i < covariates_length; ++i) {
      for(size_t j(0); j < number_types; ++j) {
        covariates_input_names[i * number_types + j] = covariates_names[i] + std::string("_") + std::to_string(j + 1);
      }
    }
  }

  Rcpp::CharacterVector log_lambda_names(Rcpp::no_init(number_types));
  for(size_t i(0); i < number_types; ++i) {
    log_lambda_names[i] = std::string("shifted_log_lambda") + std::to_string(i + 1);
  }
  Rcpp::colnames(rho_offset_rcpp) = Rcpp::wrap("rho");
  Rcpp::colnames(response_rcpp) = Rcpp::wrap("response");

  // Put all the regressors into a unique matrix that will be sent to glm / glmnet.
  Rcpp::NumericMatrix regressors;
  if(number_traits > 0) {
    regressors = Rcpp::no_init(log_lambda.nrow(),
                               log_lambda.ncol() + short_range_traits_input.ncol() + medium_range_traits_input.ncol() + short_range_joint_traits_input.ncol() + medium_range_joint_traits_input.ncol() + covariates_input.ncol());
    Rcpp::CharacterVector col_names(Rcpp::no_init(regressors.ncol()));
    for(R_xlen_t j(0); j < static_cast<R_xlen_t>(log_lambda.ncol()); ++j) {
      col_names[j] = log_lambda_names[j];
      for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
        regressors(i, j) = log_lambda(i, j);
      }
    }
    R_xlen_t index_shift(log_lambda.ncol());
    for(R_xlen_t j(0); j < static_cast<R_xlen_t>(short_range_traits_input.ncol()); ++j) {
      col_names[j + index_shift] = a_names[j];
      for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
        regressors(i, j + index_shift) = short_range_traits_input(i, j);
      }
    }
    index_shift = log_lambda.ncol() + short_range_traits_input.ncol();
    for(R_xlen_t j(0); j < static_cast<R_xlen_t>(medium_range_traits_input.ncol()); ++j) {
      col_names[j + index_shift] = b_names[j];
      for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
        regressors(i, j + index_shift) = medium_range_traits_input(i, j);
      }
    }
    index_shift = log_lambda.ncol() + short_range_traits_input.ncol() + medium_range_traits_input.ncol();
    for(R_xlen_t j(0); j < static_cast<R_xlen_t>(short_range_joint_traits_input.ncol()); ++j) {
      col_names[j + index_shift] = c_names[j];
      for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
        regressors(i, j + index_shift) = short_range_joint_traits_input(i, j);
      }
    }
    index_shift = log_lambda.ncol() + short_range_traits_input.ncol() + medium_range_traits_input.ncol() + short_range_joint_traits_input.ncol();
    for(R_xlen_t j(0); j < static_cast<R_xlen_t>(medium_range_joint_traits_input.ncol()); ++j) {
      col_names[j + index_shift] = d_names[j];
      for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
        regressors(i, j + index_shift) = medium_range_joint_traits_input(i, j);
      }
    }
    index_shift = log_lambda.ncol() + short_range_traits_input.ncol() + medium_range_traits_input.ncol() + short_range_joint_traits_input.ncol() + medium_range_joint_traits_input.ncol();
    for(R_xlen_t j(0); j < static_cast<R_xlen_t>(covariates_input.ncol()); ++j) {
      col_names[j + index_shift] = covariates_input_names[j];
      for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
        regressors(i, j + index_shift) = covariates_input(i, j);
      }
    }
    Rcpp::colnames(regressors) = col_names;
  } else {
    regressors = Rcpp::no_init(log_lambda.nrow(),
                               log_lambda.ncol() + alpha_input.ncol() + gamma_input.ncol() + covariates_input.ncol());
    Rcpp::CharacterVector col_names(Rcpp::no_init(regressors.ncol()));
    for(R_xlen_t j(0); j < static_cast<R_xlen_t>(log_lambda.ncol()); ++j) {
      col_names[j] = log_lambda_names[j];
      for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
        regressors(i, j) = log_lambda(i, j);
      }
    }
    R_xlen_t index_shift(log_lambda.ncol());
    for(R_xlen_t j(0); j < static_cast<R_xlen_t>(alpha_input.ncol()); ++j) {
      col_names[j + index_shift] = alpha_names[j];
      for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
        regressors(i, j + index_shift) = alpha_input(i, j);
      }
    }
    index_shift = log_lambda.ncol() + alpha_input.ncol();
    for(R_xlen_t j(0); j < static_cast<R_xlen_t>(gamma_input.ncol()); ++j) {
      col_names[j + index_shift] = gamma_names[j];
      for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
        regressors(i, j + index_shift) = gamma_input(i, j);
      }
    }
    index_shift = log_lambda.ncol() + alpha_input.ncol() + gamma_input.ncol();
    for(R_xlen_t j(0); j < static_cast<R_xlen_t>(covariates_input.ncol()); ++j) {
      col_names[j + index_shift] = covariates_input_names[j];
      for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
        regressors(i, j + index_shift) = covariates_input(i, j);
      }
    }
    Rcpp::colnames(regressors) = col_names;
  }


  return Rcpp::List::create(Rcpp::Named("response") = response_rcpp,
                            Rcpp::Named("offset") = rho_offset_rcpp,
                            Rcpp::Named("regressors") = regressors,
                            Rcpp::Named("shift") = shift
                            );
}

// [[Rcpp::export]]
Rcpp::List prepare_gibbsm_data(Rcpp::List configuration_list, SEXP window, Rcpp::List covariates, Rcpp::List traits, Rcpp::CharacterVector model, Rcpp::CharacterVector medium_range_model, SEXP short_range, SEXP medium_range, SEXP long_range, R_xlen_t saturation, Rcpp::NumericVector mark_range, bool approximate) {
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

  const auto cpp_window(ppjsdm::get_window_ptr_from_R_object(window, mark_range));
  // The trick below allows us to find the number of different types in the configuration.
  // That number is then used to default construct `short_range`.
  std::vector<size_t> max_points_by_type(ppjsdm::get_number_points(vector_configurations[0]));
  using size_t = typename decltype(vector_configurations)::size_type;
  for(size_t i(0); i < vector_configurations.size(); ++i) {
    const auto new_points_by_type(ppjsdm::get_number_points(vector_configurations[i]));
    for(size_t j(0); j < max_points_by_type.size(); ++j) {
      if(max_points_by_type[j] < new_points_by_type[j]) {
        max_points_by_type[j] = new_points_by_type[j];
      }
    }
  }
  // TODO: Allow for cases in which all species are not present in all configurations, i.e. number_types = max_i(max_points_by_type[i].size())
  const auto number_types(max_points_by_type.size());

  const auto sh(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(short_range, 0.1 * cpp_window->diameter(), number_types));
  if(!ppjsdm::is_symmetric_matrix(sh)) {
    Rcpp::stop("Short range interaction radius matrix is not symmetric.");
  }
  return ppjsdm::call_on_dispersion_model(model, sh, saturation, [approximate, number_types, medium_range_model, medium_range, long_range, saturation, &vector_configurations, &cpp_window, &covariates, traits, &max_points_by_type](const auto& short_papangelou) {
    const auto me(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(medium_range, 0.1 * cpp_window->diameter(), number_types));
    const auto lo(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(long_range, 0.2 * cpp_window->diameter(), number_types));
    if(!ppjsdm::is_symmetric_matrix(me) || !ppjsdm::is_symmetric_matrix(lo)) {
      Rcpp::stop("Medium or long range interaction radius matrix is not symmetric.");
    }
    return ppjsdm::call_on_medium_range_dispersion_model(medium_range_model, me, lo, saturation, [approximate, &short_papangelou, &vector_configurations, &cpp_window, &covariates, traits, &max_points_by_type](const auto& medium_papangelou) {
      if(approximate) {
        return prepare_gibbsm_data_helper<true>(vector_configurations, *cpp_window, ppjsdm::Im_list_wrapper(covariates), traits, short_papangelou, medium_papangelou, max_points_by_type);
      } else {
        return prepare_gibbsm_data_helper<false>(vector_configurations, *cpp_window, ppjsdm::Im_list_wrapper(covariates), traits, short_papangelou, medium_papangelou, max_points_by_type);
      }
    });
  });
}
