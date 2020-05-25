#include <RcppArmadillo.h>
#include <Rinternals.h>

#include "configuration/configuration_manipulation.hpp"
#include "configuration/configuration_wrapper.hpp"
#include "configuration/get_number_points.hpp"

#include "point/point_manipulation.hpp"

#include "saturated_model/exponential_family_model.hpp"
#include "saturated_model/saturated_model.hpp"

#include "simulation/rbinomialpp_single.hpp"

#include "utility/construct_if_missing.hpp"
#include "utility/im_wrapper.hpp"
#include "utility/is_symmetric_matrix.hpp"
#include "utility/lightweight_matrix.hpp"
#include "utility/window.hpp"

#include <algorithm> // std::remove_if
#include <cmath> // std::log
#include <iterator> // std::next
#include <string> // std::string, std::to_string
#include <vector> // std::vector

#ifdef _OPENMP
#include <omp.h>
#endif

template<typename Configuration, typename Model>
Rcpp::NumericMatrix compute_vcov_helper(const Configuration& configuration,
                                         const ppjsdm::Im_list_wrapper& covariates,
                                         const ppjsdm::Saturated_model& dispersion_model,
                                         const ppjsdm::Saturated_model& medium_dispersion_model,
                                         const Model& model,
                                         int number_types,
                                         double rho,
                                         Rcpp::NumericVector coefficients_vector,
                                         Rcpp::NumericMatrix regressors,
                                         Rcpp::List data_list,
                                         bool estimate_alpha,
                                         bool estimate_gamma) {
  using size_t = ppjsdm::size_t<Configuration>;

  const size_t total_points(ppjsdm::size(configuration));
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
  const size_t total_parameters(index_start_covariates + number_types * covariates.size());

  const auto x(Rcpp::as<Rcpp::NumericVector>(data_list["x"]));
  const auto y(Rcpp::as<Rcpp::NumericVector>(data_list["y"]));
  const auto type(Rcpp::as<Rcpp::IntegerVector>(data_list["type"]));
  const auto mark(Rcpp::as<Rcpp::NumericVector>(data_list["mark"]));

  arma::mat A(total_parameters, total_parameters);
  arma::mat S(total_parameters, total_parameters);
  for(size_t k1(0); k1 < total_parameters; ++k1) {
    for(size_t k2(0); k2 < total_parameters; ++k2) {
      A(k1, k2) = 0.;
      S(k1, k2) = 0.;
    }
  }

  std::vector<double> papangelou(regressors.nrow());
  double kappa(0);

  // The lines below set A to G2 from Baddeley et al
  for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
    double inner_product(0);
    for(R_xlen_t j(0); j < coefficients_vector.size(); ++j) {
      inner_product += coefficients_vector[j] * regressors(i, j);
    }
    const auto papangelou_value(std::exp(inner_product));
    papangelou[i] = papangelou_value;
    kappa += 1. / (papangelou_value + rho);

    const auto constant(papangelou_value / ((papangelou_value + rho) * (papangelou_value + rho)));
    for(size_t k1(0); k1 < total_parameters; ++k1) {
      const auto value(regressors(i, k1) * constant);
      for(size_t k2(0); k2 < total_parameters; ++k2) {
        A(k1, k2) += value;
        S(k1, k2) += value * regressors(i, k2);
      }
    }
  }
  A *= A.t();
  A *= -1. / kappa;
  for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
    const auto constant(papangelou[i] * papangelou[i] / ((papangelou[i] + rho) * (papangelou[i] + rho) * (papangelou[i] + rho)));
    for(size_t k1(0); k1 < total_parameters; ++k1) {
      A(k1, k1) += regressors(i, k1) * regressors(i, k1) * constant;
      for(size_t k2(k1 + 1); k2 < total_parameters; ++k2) {
        A(k1, k2) += regressors(i, k1) * regressors(i, k2) * constant;
      }
    }
  }
  A /= rho;

  // At this point, A is equal to G2. Next, add A1.
  for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
    const auto constant(papangelou[i] / ((papangelou[i] + rho) * (papangelou[i] + rho) * (papangelou[i] + rho)));
    for(size_t k1(0); k1 < total_parameters; ++k1) {
      A(k1, k1) += regressors(i, k1) * regressors(i, k1) * constant;
      for(size_t k2(k1 + 1); k2 < total_parameters; ++k2) {
        A(k1, k2) += regressors(i, k1) * regressors(i, k2) * constant;
        A(k2, k1) = A(k1, k2);
      }
    }
  }

  // Finally, add A2 and A3.
  for(size_t i(0); i < total_points; ++i) {
    const int type_i(ppjsdm::get_type(configuration[i]));

    std::vector<double> cov_i(covariates.size());
    for(R_xlen_t k(0); k < covariates.size(); ++k) {
      cov_i[k] = covariates[k](configuration[i]);
    }

    std::vector<double> t_over_papangelou_i(total_parameters);
    for(R_xlen_t l(0); l < regressors.nrow(); ++l) {
      if(ppjsdm::is_equal(configuration[i], ppjsdm::Marked_point(x[l], y[l], type[l] - 1, mark[l]))) {
        const auto one_over_papangelou_plus_rho(1. / (papangelou[l] + rho));
        for(decltype(t_over_papangelou_i.size()) fill(0); fill < t_over_papangelou_i.size(); ++fill) {
            t_over_papangelou_i[fill] = regressors(l, fill) * one_over_papangelou_plus_rho;
        }
        break;
      }
    }
    // if(t_over_papangelou_i.size() == 0) {
    //   Rcpp::stop("Did not find the current point in t_over_papangelou");
    // }

    Configuration configuration_without_i(configuration);
    ppjsdm::remove_point(configuration_without_i, configuration[i]);
    const auto papangelou_i_minus_i(model.compute_papangelou(configuration[i], configuration_without_i));

    for(size_t j(i + 1); j < total_points; ++j) {
      const int type_j(ppjsdm::get_type(configuration[j]));
      std::vector<double> t_over_papangelou_j(total_parameters);
      for(R_xlen_t l(0); l < regressors.nrow(); ++l) {
        if(ppjsdm::is_equal(configuration[j], ppjsdm::Marked_point(x[l], y[l], type[l] - 1, mark[l]))) {
          const auto one_over_papangelou_plus_rho(1. / (papangelou[l] + rho));
          for(decltype(t_over_papangelou_j.size()) fill(0); fill < t_over_papangelou_j.size(); ++fill) {
            t_over_papangelou_j[fill] = regressors(l, fill) * one_over_papangelou_plus_rho;
          }
          break;
        }
      }

      Configuration configuration_without_ij(configuration_without_i);
      ppjsdm::remove_point(configuration_without_ij, configuration[j]);

      const auto papangelou_i_minus_two(model.compute_papangelou(configuration[i], configuration_without_ij));
      const auto papangelou_j_minus_two(model.compute_papangelou(configuration[j], configuration_without_ij));

      // TODO: Don't need to compute all of these depending on estimate_alpha / estimate_gamma
      const auto short_i(ppjsdm::compute_dispersion(dispersion_model, configuration[i], number_types, configuration_without_ij));
      const auto medium_i(ppjsdm::compute_dispersion(medium_dispersion_model, configuration[i], number_types, configuration_without_ij));

      const auto short_j(ppjsdm::compute_dispersion(dispersion_model, configuration[j], number_types, configuration_without_ij));
      const auto medium_j(ppjsdm::compute_dispersion(medium_dispersion_model, configuration[j], number_types, configuration_without_ij));

      // TODO: Reuse regressors?
      std::vector<double> cov_j(covariates.size());
      for(R_xlen_t k(0); k < covariates.size(); ++k) {
        cov_j[k] = covariates[k](configuration[j]);
      }

      std::vector<double> t_i(total_parameters);
      std::vector<double> t_j(total_parameters);

      size_t current_index(0);
      for(int k1(0); k1 < number_types; ++k1) {
        if(k1 == type_i) {
          t_i[k1] = 1.;

          for(int k2(0); k2 < covariates.size(); ++k2) {
            t_i[index_start_covariates + k2 * number_types + k1] = cov_i[k2];
          }

          size_t filling(current_index);
          for(int k2(k1); k2 < number_types; ++k2) {
            if(estimate_alpha) {
              t_i[number_types + filling] = short_i[k2];
            }
            if(estimate_gamma) {
              t_i[index_start_gamma + filling] = medium_i[k2];
            }
            ++filling;
          }
        } else {
          size_t filling(current_index);
          for(int k2(k1); k2 < number_types; ++k2) {
            if(k2 == type_i) {
              if(estimate_alpha) {
                t_i[number_types + filling] = short_i[k1];
              }
              if(estimate_gamma) {
                t_i[index_start_gamma + filling] = medium_i[k1];
              }
            }
            ++filling;
          }
        }

        if(k1 == type_j) {
          t_j[k1] = 1.;

          for(int k2(0); k2 < covariates.size(); ++k2) {
            t_j[index_start_covariates + k2 * number_types + k1] = cov_j[k2];
          }

          for(int k2(k1); k2 < number_types; ++k2) {
            if(estimate_alpha) {
              t_j[number_types + current_index] = short_j[k2];
            }
            if(estimate_gamma) {
              t_j[index_start_gamma + current_index] = medium_j[k2];
            }
            ++current_index;
          }
        } else {
          for(int k2(k1); k2 < number_types; ++k2) {
            if(k2 == type_j) {
              if(estimate_alpha) {
                t_j[number_types + current_index] = short_j[k1];
              }
              if(estimate_gamma) {
                t_j[index_start_gamma + current_index] = medium_j[k1];
              }
            }
            ++current_index;
          }
        }
      }

      const auto constant(2. * (papangelou_i_minus_two / papangelou_i_minus_i - 1.) / ((papangelou_i_minus_two + rho) * (papangelou_j_minus_two + rho)));
      for(size_t k1(0); k1 < total_parameters; ++k1) {
        for(size_t k2(0); k2 < total_parameters; ++k2) {
          const auto A2_summand(t_i[k1] * t_j[k2] * constant);
          const auto A3_summand(2. * (t_over_papangelou_i[k1] - t_i[k1] / (papangelou_i_minus_two + rho)) * (t_over_papangelou_j[k2] - t_j[k2] / (papangelou_j_minus_two + rho)));
          A(k1, k2) += A2_summand + A3_summand;
        }
      }
    }
  }

  const arma::mat S_inv(arma::inv(S));
  A = S_inv * A * S_inv;

  // Set names
  Rcpp::CharacterVector col_names(Rcpp::no_init(A.n_cols));
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

  Rcpp::NumericMatrix vcov(Rcpp::wrap(A));
  Rcpp::colnames(vcov) = col_names;
  Rcpp::rownames(vcov) = col_names;

  return vcov;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_vcov(SEXP configuration,
                                 Rcpp::List covariates,
                                 Rcpp::CharacterVector model,
                                 Rcpp::CharacterVector medium_range_model,
                                 SEXP short_range,
                                 SEXP medium_range,
                                 SEXP long_range,
                                 R_xlen_t saturation,
                                 Rcpp::NumericMatrix alpha,
                                 Rcpp::NumericVector beta0,
                                 Rcpp::NumericMatrix beta,
                                 Rcpp::NumericMatrix gamma,
                                 double rho,
                                 Rcpp::NumericVector coefficients_vector,
                                 Rcpp::NumericMatrix regressors,
                                 Rcpp::List data_list,
                                 bool estimate_alpha,
                                 bool estimate_gamma) {
  // Construct std::vector of configurations.
  const ppjsdm::Configuration_wrapper wrapped_configuration(Rcpp::wrap(configuration));
  const auto length_configuration(ppjsdm::size(wrapped_configuration));

  // Convert configuration to std::vector in order for parallelised version to work.
  std::vector<ppjsdm::Marked_point> vector_configuration(length_configuration);
  for(decltype(ppjsdm::size(wrapped_configuration)) j(0); j < length_configuration; ++j) {
    vector_configuration[j] = wrapped_configuration[j];
  }

  const auto number_types(beta0.size());
  const auto dispersion(ppjsdm::Saturated_model(model, short_range, saturation));
  const auto medium_range_dispersion(ppjsdm::Saturated_model(medium_range_model, medium_range, long_range, saturation));

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

  return compute_vcov_helper(vector_configuration,
                              ppjsdm::Im_list_wrapper(covariates),
                              dispersion,
                              medium_range_dispersion,
                              exponential_model,
                              number_types,
                              rho,
                              coefficients_vector,
                              regressors,
                              data_list,
                              estimate_alpha,
                              estimate_gamma);
}
