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
#include "utility/gibbsm_helper_functions.hpp"
#include "utility/im_wrapper.hpp"
#include "utility/is_symmetric_matrix.hpp"
#include "utility/lightweight_matrix.hpp"
#include "utility/window.hpp"

#include <algorithm> // std::remove_if
#include <cmath> // std::log
#include <iterator> // std::next
#include <string> // std::string, std::to_string
#include <vector> // std::vector

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

  const auto number_parameters_struct(ppjsdm::get_number_parameters(number_types, covariates.size(), estimate_alpha, estimate_gamma));
  const auto index_start_gamma(number_parameters_struct.index_start_gamma);
  const auto index_start_covariates(number_parameters_struct.index_start_covariates);
  const auto total_parameters(number_parameters_struct.total_parameters);

  const auto response(Rcpp::as<Rcpp::IntegerVector>(data_list["response"]));
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

  // Fill linear values
  std::vector<std::vector<double>> covariates_on_configuration(regressors.nrow());
  std::vector<std::vector<double>> t_over_papangelou_on_configuration(regressors.nrow());
  for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
    if(response[i] == 1) {
      const ppjsdm::Marked_point point_i(x[i], y[i], type[i] - 1, mark[i]);

      std::vector<double> cov(covariates.size());
      for(R_xlen_t k(0); k < covariates.size(); ++k) {
        cov[k] = regressors(i, index_start_covariates + k * number_types + ppjsdm::get_type(point_i));
      }
      covariates_on_configuration[i] = cov;

      std::vector<double> t_over_papangelou(total_parameters);
      const auto one_over_papangelou_plus_rho(1. / (papangelou[i] + rho));
      for(decltype(t_over_papangelou.size()) fill(0); fill < t_over_papangelou.size(); ++fill) {
        t_over_papangelou[fill] = regressors(i, fill) * one_over_papangelou_plus_rho;
      }
      t_over_papangelou_on_configuration[i] = t_over_papangelou;
    }
  }

  // Finally, add A2 and A3.
  for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
    if(response[i] == 1) {
      const ppjsdm::Marked_point point_i(x[i], y[i], type[i] - 1, mark[i]);

      Configuration configuration_without_i(configuration);
      ppjsdm::remove_point(configuration_without_i, point_i);
      const auto papangelou_i_minus_i(model.compute_papangelou(point_i, configuration_without_i));

      for(R_xlen_t j(i + 1); j < regressors.nrow(); ++j) {
        if(response[j] == 1) {
          const ppjsdm::Marked_point point_j(x[j], y[j], type[j] - 1, mark[j]);

          Configuration configuration_without_ij(configuration_without_i);
          ppjsdm::remove_point(configuration_without_ij, point_j);

          const auto papangelou_i_minus_two(model.compute_papangelou(point_i, configuration_without_ij));
          const auto papangelou_j_minus_two(model.compute_papangelou(point_j, configuration_without_ij));

          using dispersion_t = decltype(ppjsdm::compute_dispersion(dispersion_model, point_i, number_types, configuration_without_ij));
          dispersion_t short_i, short_j, medium_i, medium_j;
          if(estimate_alpha) {
            short_i = ppjsdm::compute_dispersion(dispersion_model, point_i, number_types, configuration_without_ij);
            short_j = ppjsdm::compute_dispersion(dispersion_model, point_j, number_types, configuration_without_ij);
          }
          if(estimate_gamma) {
            medium_i = ppjsdm::compute_dispersion(medium_dispersion_model, point_i, number_types, configuration_without_ij);
            medium_j = ppjsdm::compute_dispersion(medium_dispersion_model, point_j, number_types, configuration_without_ij);
          }

          std::vector<double> t_i(total_parameters);
          std::vector<double> t_j(total_parameters);

          size_t current_index(0);
          for(int k1(0); k1 < number_types; ++k1) {
            if(k1 == ppjsdm::get_type(point_i)) {
              t_i[k1] = 1.;

              for(int k2(0); k2 < covariates.size(); ++k2) {
                t_i[index_start_covariates + k2 * number_types + k1] = covariates_on_configuration[i][k2];
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
                if(k2 == ppjsdm::get_type(point_i)) {
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

            if(k1 == ppjsdm::get_type(point_j)) {
              t_j[k1] = 1.;

              for(int k2(0); k2 < covariates.size(); ++k2) {
                t_j[index_start_covariates + k2 * number_types + k1] = covariates_on_configuration[j][k2];
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
                if(k2 == ppjsdm::get_type(point_j)) {
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
              const auto A3_summand(2. * (t_over_papangelou_on_configuration[i][k1] - t_i[k1] / (papangelou_i_minus_two + rho)) * (t_over_papangelou_on_configuration[j][k2] - t_j[k2] / (papangelou_j_minus_two + rho)));
              A(k1, k2) += A2_summand + A3_summand;
            }
          }
        }
      }
    }
  }

  const arma::mat S_inv(arma::inv(S));
  A = S_inv * A * S_inv;

  const auto col_names(make_model_coloumn_names(covariates, number_types, estimate_alpha, estimate_gamma));

  Rcpp::NumericMatrix rcpp_matrix(Rcpp::wrap(A));
  Rcpp::colnames(rcpp_matrix) = col_names;
  Rcpp::rownames(rcpp_matrix) = col_names;

  return rcpp_matrix;
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
