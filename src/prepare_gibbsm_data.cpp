#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_manipulation.hpp"
#include "configuration/configuration_wrapper.hpp"
#include "configuration/get_number_points.hpp"

#include "point/point_manipulation.hpp"

#include "saturated_model/saturated_model.hpp"

#include "simulation/rppp_single.hpp"

#include "utility/construct_if_missing.hpp"
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

inline auto get_number_parameters(int number_types,
                                  unsigned long long int covariates_size,
                                  bool estimate_alpha,
                                  bool estimate_gamma) {
  using size_t = unsigned long long int;
  struct Number_parameters_struct {
    size_t index_start_gamma;
    size_t index_start_covariates;
    size_t total_parameters;
  };

  size_t index_start_gamma(0);
  size_t index_start_covariates(0);
  if(estimate_alpha) {
    if(estimate_gamma) {
      index_start_gamma = number_types + number_types * (number_types + 1) / 2;
      index_start_covariates = number_types * (2 + number_types);
    } else {
      index_start_covariates = number_types + number_types * (number_types + 1) / 2;
    }
  } else {
    if(estimate_gamma) {
      index_start_gamma = number_types;
      index_start_covariates = number_types + number_types * (number_types + 1) / 2;
    } else {
      index_start_covariates = number_types;
    }
  }
  const size_t total_parameters(index_start_covariates + number_types * covariates_size);
  return Number_parameters_struct{index_start_gamma, index_start_covariates, total_parameters};
}

inline auto make_model_coloumn_names(const ppjsdm::Im_list_wrapper& covariates,
                                     int number_types,
                                     bool estimate_alpha,
                                     bool estimate_gamma) {
  Rcpp::CharacterVector col_names(Rcpp::no_init(get_number_parameters(number_types, covariates.size(), estimate_alpha, estimate_gamma).total_parameters));
  for(R_xlen_t j(0); j < number_types; ++j) {
    col_names[j] = std::string("log_lambda") + std::to_string(j + 1);
  }
  R_xlen_t index_shift(number_types);
  if(estimate_alpha) {
    for(R_xlen_t k1(0); k1 < number_types; ++k1) {
      for(R_xlen_t k2(k1); k2 < number_types; ++k2) {
        col_names[index_shift] = std::string("alpha_") + std::to_string(k1 + 1) + std::string("_") + std::to_string(k2 + 1);
        ++index_shift;
      }
    }
  }
  if(estimate_gamma) {
    for(R_xlen_t k1(0); k1 < number_types; ++k1) {
      for(R_xlen_t k2(k1); k2 < number_types; ++k2) {
        col_names[index_shift] = std::string("gamma_") + std::to_string(k1 + 1) + std::string("_") + std::to_string(k2 + 1);
        ++index_shift;
      }
    }
  }
  for(R_xlen_t j(0); j < static_cast<R_xlen_t>(covariates.size()); ++j) {
    for(R_xlen_t k(0); k < number_types; ++k) {
      col_names[index_shift] = covariates.names()[j] + std::string("_") + std::to_string(k + 1);
      ++index_shift;
    }
  }

  return col_names;
}

template<bool Approximate, typename Configuration>
Rcpp::List prepare_gibbsm_data_helper(const std::vector<Configuration>& configuration_list,
                                      const ppjsdm::Window& window,
                                      const ppjsdm::Im_list_wrapper& covariates,
                                      Rcpp::List traits,
                                      const ppjsdm::Saturated_model& dispersion_model,
                                      const ppjsdm::Saturated_model& medium_dispersion_model,
                                      unsigned long long int max_points_in_any_type,
                                      R_xlen_t ndummy,
                                      bool estimate_alpha,
                                      bool estimate_gamma,
                                      int number_types) {
  using size_t = ppjsdm::size_t<Configuration>;

  // Sample the dummy points D.
  // This choice of rho is the guideline from the Baddeley et al. paper, see p. 8 therein.
  std::vector<double> rho_times_volume;
  size_t length_D(0);
  if(ndummy == 0) {
    rho_times_volume = std::vector<double>(number_types);
    const auto four_times_max_points_in_any_type(4 * max_points_in_any_type);
    if(four_times_max_points_in_any_type < 500) {
      std::fill(rho_times_volume.begin(), rho_times_volume.end(), four_times_max_points_in_any_type);
      length_D = number_types * four_times_max_points_in_any_type;
    } else {
      std::fill(rho_times_volume.begin(), rho_times_volume.end(), 500);
      length_D = number_types * 500;
    }
  } else {
    rho_times_volume = std::vector<double>(number_types);
    std::fill(rho_times_volume.begin(), rho_times_volume.end(), ndummy);
    length_D = number_types * ndummy;
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
      const auto d(ppjsdm::compute_dispersion<Approximate>(dispersion_model, D[i], number_types, configuration_list[j]));
      const auto e(ppjsdm::compute_dispersion<Approximate>(medium_dispersion_model, D[i], number_types, configuration_list[j]));
      results_private.emplace_back(false, ppjsdm::get_type(D[i]), std::move(d), std::move(e), cov, ppjsdm::get_x(D[i]), ppjsdm::get_y(D[i]), ppjsdm::get_mark(D[i]));
    }
  }
#pragma omp critical
  precomputed_results.insert(precomputed_results.end(), results_private.begin(), results_private.end());
  }

  const auto total_points(precomputed_results.size());
  //const size_t number_traits(traits.size());

  const auto number_parameters_struct(get_number_parameters(number_types, covariates.size(), estimate_alpha, estimate_gamma));
  const auto index_start_gamma(number_parameters_struct.index_start_gamma);
  const auto index_start_covariates(number_parameters_struct.index_start_covariates);
  const auto total_parameters(number_parameters_struct.total_parameters);

  // Default-initialise the data
  std::vector<int> response(total_points);
  std::vector<double> x(total_points);
  std::vector<double> y(total_points);
  std::vector<double> mark(total_points);
  std::vector<double> type(total_points);
  std::vector<double> rho_offset(total_points);
  Rcpp::NumericMatrix regressors(total_points, total_parameters);
  //ppjsdm::Lightweight_matrix<double> log_lambda(total_points, number_types);
  //ppjsdm::Lightweight_matrix<double> alpha_input(total_points, number_types * (number_types + 1) / 2);
  //ppjsdm::Lightweight_matrix<double> gamma_input(total_points, number_types * (number_types + 1) / 2);
  //ppjsdm::Lightweight_matrix<double> covariates_input(total_points, covariates_length * number_types);
  //ppjsdm::Lightweight_matrix<double> short_range_traits_input(total_points, 1 + number_traits);
  //ppjsdm::Lightweight_matrix<double> short_range_joint_traits_input(total_points, 1 + number_traits);
  //ppjsdm::Lightweight_matrix<double> medium_range_traits_input(total_points, 1 + number_traits);
  //ppjsdm::Lightweight_matrix<double> medium_range_joint_traits_input(total_points, 1 + number_traits);


  // Fill the regressors, response, offset and shift with what we precomputed.
  for(size_t i(0); i < total_points; ++i) {
    if(precomputed_results[i].is_in_configuration) {
      response[i] = 1;
    }
    x[i] = precomputed_results[i].x;
    y[i] = precomputed_results[i].y;
    mark[i] = precomputed_results[i].mark;

    const int type_index(precomputed_results[i].type);
    type[i] = type_index;

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

    size_t index(0);
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
          if(estimate_alpha) {
            regressors(i, number_types + index) = precomputed_results[i].dispersion[k];
          }
          if(estimate_gamma) {
            regressors(i, index_start_gamma + index) = precomputed_results[i].medium_dispersion[k];
          }
          ++index;
        }
      } else {
        // fill alpha & gamma
        for(int k(j); k < number_types; ++k) {
          if(k == type_index) {
            if(estimate_alpha) {
              regressors(i, number_types + index) = precomputed_results[i].dispersion[j];
            }
            if(estimate_gamma) {
              regressors(i, index_start_gamma + index) = precomputed_results[i].medium_dispersion[j];
            }
          }
          ++index;
        }
      }
    }
  }

  // Convert response and rho_offset to Rcpp objects.
  Rcpp::IntegerMatrix response_rcpp(Rcpp::no_init(response.size(), 1));
  Rcpp::NumericMatrix x_rcpp(Rcpp::no_init(response.size(), 1));
  Rcpp::NumericMatrix y_rcpp(Rcpp::no_init(response.size(), 1));
  Rcpp::NumericMatrix mark_rcpp(Rcpp::no_init(response.size(), 1));
  Rcpp::IntegerMatrix type_rcpp(Rcpp::no_init(response.size(), 1));
  Rcpp::NumericMatrix rho_offset_rcpp(Rcpp::no_init(response.size(), 1));
  for(R_xlen_t i(0); i < static_cast<R_xlen_t>(response.size()); ++i) {
    response_rcpp(i, 0) = response[i];
    x_rcpp(i, 0) = x[i];
    y_rcpp(i, 0) = y[i];
    mark_rcpp(i, 0) = mark[i];
    type_rcpp(i, 0) = type[i] + 1;
    rho_offset_rcpp(i, 0) = rho_offset[i];
  }

  // Set names.
  // Rcpp::CharacterVector alpha_names(Rcpp::no_init(alpha_input.ncol()));
  // Rcpp::CharacterVector gamma_names(Rcpp::no_init(alpha_input.ncol()));
  //
  // size_t index(0);
  // for(int j(0); j < number_types; ++j) {
  //   for(int k(j); k < number_types; ++k) {
  //     alpha_names[index] = std::string("alpha_") + std::to_string(j + 1) + std::string("_") + std::to_string(k + 1);
  //     gamma_names[index++] = std::string("gamma_") + std::to_string(j + 1) + std::string("_") + std::to_string(k + 1);
  //   }
  // }

  // Rcpp::CharacterVector short_range_direct_names(Rcpp::no_init(short_range_traits_input.ncol()));
  // Rcpp::CharacterVector medium_range_direct_names(Rcpp::no_init(medium_range_traits_input.ncol()));
  // Rcpp::CharacterVector short_range_joint_names(Rcpp::no_init(short_range_joint_traits_input.ncol()));
  // Rcpp::CharacterVector medium_range_joint_names(Rcpp::no_init(medium_range_joint_traits_input.ncol()));
  // for(size_t i(0); i < short_range_traits_input.ncol(); ++i) {
  //   short_range_direct_names[i] = std::string("short_range_direct_") + (i == 0 ? std::string("intercept") : std::to_string(i));
  //   medium_range_direct_names[i] = std::string("medium_range_direct_") + (i == 0 ? std::string("intercept") : std::to_string(i));
  //   short_range_joint_names[i] = std::string("short_range_joint_") + (i == 0 ? std::string("intercept") : std::to_string(i));
  //   medium_range_joint_names[i] = std::string("medium_range_joint_") + (i == 0 ? std::string("intercept") : std::to_string(i));
  // }

  // Rcpp::CharacterVector covariates_input_names(Rcpp::no_init(covariates_length * number_types));
  // if(covariates_length > 0) {
  //   const auto covariates_names(covariates.names());
  //   for(size_t i(0); i < covariates_length; ++i) {
  //     for(int j(0); j < number_types; ++j) {
  //       covariates_input_names[i * number_types + j] = covariates_names[i] + std::string("_") + std::to_string(j + 1);
  //     }
  //   }
  // }
  //
  // Rcpp::CharacterVector log_lambda_names(Rcpp::no_init(number_types));
  // for(int i(0); i < number_types; ++i) {
  //   log_lambda_names[i] = std::string("log_lambda") + std::to_string(i + 1);
  // }
  Rcpp::colnames(rho_offset_rcpp) = Rcpp::wrap("rho");
  Rcpp::colnames(response_rcpp) = Rcpp::wrap("response");
  Rcpp::colnames(x_rcpp) = Rcpp::wrap("x");
  Rcpp::colnames(y_rcpp) = Rcpp::wrap("y");
  Rcpp::colnames(mark_rcpp) = Rcpp::wrap("mark");
  Rcpp::colnames(type_rcpp) = Rcpp::wrap("type");

  // if(number_traits > 0) {
  //   regressors = Rcpp::no_init(log_lambda.nrow(),
  //                              log_lambda.ncol() + short_range_traits_input.ncol() + medium_range_traits_input.ncol() + short_range_joint_traits_input.ncol() + medium_range_joint_traits_input.ncol() + covariates_input.ncol());
  //   Rcpp::CharacterVector col_names(Rcpp::no_init(regressors.ncol()));
  //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(log_lambda.ncol()); ++j) {
  //     col_names[j] = log_lambda_names[j];
  //     for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
  //       regressors(i, j) = log_lambda(i, j);
  //     }
  //   }
  //   R_xlen_t index_shift(log_lambda.ncol());
  //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(short_range_traits_input.ncol()); ++j) {
  //     col_names[j + index_shift] = short_range_direct_names[j];
  //     for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
  //       regressors(i, j + index_shift) = short_range_traits_input(i, j);
  //     }
  //   }
  //   index_shift = log_lambda.ncol() + short_range_traits_input.ncol();
  //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(medium_range_traits_input.ncol()); ++j) {
  //     col_names[j + index_shift] = medium_range_direct_names[j];
  //     for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
  //       regressors(i, j + index_shift) = medium_range_traits_input(i, j);
  //     }
  //   }
  //   index_shift = log_lambda.ncol() + short_range_traits_input.ncol() + medium_range_traits_input.ncol();
  //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(short_range_joint_traits_input.ncol()); ++j) {
  //     col_names[j + index_shift] = short_range_joint_names[j];
  //     for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
  //       regressors(i, j + index_shift) = short_range_joint_traits_input(i, j);
  //     }
  //   }
  //   index_shift = log_lambda.ncol() + short_range_traits_input.ncol() + medium_range_traits_input.ncol() + short_range_joint_traits_input.ncol();
  //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(medium_range_joint_traits_input.ncol()); ++j) {
  //     col_names[j + index_shift] = medium_range_joint_names[j];
  //     for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
  //       regressors(i, j + index_shift) = medium_range_joint_traits_input(i, j);
  //     }
  //   }
  //   index_shift = log_lambda.ncol() + short_range_traits_input.ncol() + medium_range_traits_input.ncol() + short_range_joint_traits_input.ncol() + medium_range_joint_traits_input.ncol();
  //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(covariates_input.ncol()); ++j) {
  //     col_names[j + index_shift] = covariates_input_names[j];
  //     for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
  //       regressors(i, j + index_shift) = covariates_input(i, j);
  //     }
  //   }
  //   Rcpp::colnames(regressors) = col_names;
  // } else {
    // Rcpp::CharacterVector col_names(Rcpp::no_init(regressors.ncol()));
    // for(R_xlen_t j(0); j < static_cast<R_xlen_t>(log_lambda.ncol()); ++j) {
    //   col_names[j] = log_lambda_names[j];
    // }
    // R_xlen_t index_shift(log_lambda.ncol());
    // if(estimate_alpha) {
    //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(alpha_input.ncol()); ++j) {
    //     col_names[j + index_shift] = alpha_names[j];
    //   }
    //   index_shift += alpha_input.ncol();
    // }
    // if(estimate_gamma) {
    //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(gamma_input.ncol()); ++j) {
    //     col_names[j + index_shift] = gamma_names[j];
    //   }
    //   index_shift += gamma_input.ncol();
    // }
    // for(R_xlen_t j(0); j < static_cast<R_xlen_t>(covariates_input.ncol()); ++j) {
    //   col_names[j + index_shift] = covariates_input_names[j];
    // }
    const auto col_names(make_model_coloumn_names(covariates, number_types, estimate_alpha, estimate_gamma));
    Rcpp::colnames(regressors) = col_names;
 // }


  return Rcpp::List::create(Rcpp::Named("response") = response_rcpp,
                            Rcpp::Named("x") = x_rcpp,
                            Rcpp::Named("y") = y_rcpp,
                            Rcpp::Named("mark") = mark_rcpp,
                            Rcpp::Named("type") = type_rcpp,
                            Rcpp::Named("offset") = rho_offset_rcpp,
                            Rcpp::Named("regressors") = regressors,
                            Rcpp::Named("shift") = shift
                            );
}

// [[Rcpp::export]]
Rcpp::List prepare_gibbsm_data(Rcpp::List configuration_list,
                               SEXP window,
                               Rcpp::List covariates,
                               Rcpp::List traits,
                               Rcpp::CharacterVector model,
                               Rcpp::CharacterVector medium_range_model,
                               SEXP short_range,
                               SEXP medium_range,
                               SEXP long_range,
                               R_xlen_t saturation,
                               Rcpp::NumericVector mark_range,
                               bool approximate,
                               R_xlen_t ndummy,
                               bool estimate_alpha,
                               bool estimate_gamma) {
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
  // The trick below allows us to find the number of different types in the configuration.
  // That number is then used to default construct `short_range`.
  const auto max_points_by_type(ppjsdm::get_number_points(vector_configurations[0]));
  const auto number_types(max_points_by_type.size());
  auto max_points_in_any_type(*std::max_element(max_points_by_type.begin(), max_points_by_type.end()));

  using size_t = typename decltype(vector_configurations)::size_type;
  for(size_t i(0); i < vector_configurations.size(); ++i) {
    const auto new_points_by_type(ppjsdm::get_number_points(vector_configurations[i]));
    for(size_t j(0); j < max_points_by_type.size(); ++j) {
      if(max_points_in_any_type < new_points_by_type[j]) {
        max_points_in_any_type = new_points_by_type[j];
      }
    }
  }
  const auto dispersion(ppjsdm::Saturated_model(model, short_range, saturation));
  const auto medium_range_dispersion(ppjsdm::Saturated_model(medium_range_model, medium_range, long_range, saturation));
  if(approximate) {
    return prepare_gibbsm_data_helper<true>(vector_configurations, cpp_window, ppjsdm::Im_list_wrapper(covariates), traits, dispersion, medium_range_dispersion, max_points_in_any_type, ndummy, estimate_alpha, estimate_gamma, number_types);
  } else {
    return prepare_gibbsm_data_helper<false>(vector_configurations, cpp_window, ppjsdm::Im_list_wrapper(covariates), traits, dispersion, medium_range_dispersion, max_points_in_any_type, ndummy, estimate_alpha, estimate_gamma, number_types);
  }
}
