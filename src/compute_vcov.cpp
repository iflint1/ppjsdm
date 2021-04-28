#include <RcppArmadillo.h>
#include <Rinternals.h>

#include "configuration/configuration_manipulation.hpp"
#include "configuration/configuration_wrapper.hpp"
#include "configuration/get_number_points.hpp"

#include "point/point_manipulation.hpp"

#include "saturated_model/compute_dispersion_vcov.hpp"
#include "saturated_model/exponential_family_model.hpp"
#include "saturated_model/saturated_model.hpp"

#include "simulation/rbinomialpp_single.hpp"

#include "utility/construct_if_missing.hpp"
#include "utility/flatten_strict_upper_triangular.hpp"
#include "utility/gibbsm_helper_functions.hpp"
#include "utility/im_wrapper.hpp"
#include "utility/is_symmetric_matrix.hpp"
#include "utility/lightweight_matrix.hpp"
#include "utility/window.hpp"

#include <algorithm> // std::remove_if
#include <cmath> // std::floor, std::log, std::sqrt
#include <iterator> // std::next
#include <string> // std::string, std::to_string
#include <type_traits> // std::remove_cv_t
#include <vector> // std::vector

template<typename Configuration>
Rcpp::NumericMatrix compute_vcov_helper(const Configuration& configuration,
                                        const ppjsdm::Im_list_wrapper& covariates,
                                        const ppjsdm::Saturated_model& dispersion_model,
                                        const ppjsdm::Saturated_model& medium_dispersion_model,
                                        int number_types,
                                        double rho,
                                        Rcpp::NumericVector coefficients_vector,
                                        Rcpp::NumericMatrix regressors,
                                        Rcpp::List data_list,
                                        Rcpp::LogicalMatrix estimate_alpha,
                                        Rcpp::LogicalMatrix estimate_gamma,
                                        int nthreads) {

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

  // Do you need to compute some of the alphas/gammas?
  bool compute_some_alphas(false);
  bool compute_some_gammas(false);
  for(int i(0); i < number_types; ++i) {
    for(int j(i); j < number_types; ++j) {
      if(estimate_alpha(i, j)) {
        compute_some_alphas = true;
      }
      if(estimate_gamma(i, j)) {
        compute_some_gammas = true;
      }
    }
  }

  // Parallel computations done before adding A2 and A3
  const auto precomputed_size(regressors.nrow() * (regressors.nrow() - 1) / 2);

  std::vector<R_xlen_t> non_zero_ij_responses;
  for(std::remove_cv_t<decltype(precomputed_size)> k(0); k < precomputed_size; ++k) {
    const auto pr(ppjsdm::decode_linear(k, regressors.nrow()));
    if(response[pr.first] == 1 && response[pr.second] == 1) {
      non_zero_ij_responses.emplace_back(k);
    }
  }

  using precomputation_t = decltype(ppjsdm::compute_dispersion_for_vcov(dispersion_model, number_types, configuration, nthreads));
  precomputation_t short_computation, medium_computation;
  if(compute_some_alphas) {
    short_computation = ppjsdm::compute_dispersion_for_vcov(dispersion_model, number_types, configuration, nthreads);
  }
  if(compute_some_gammas) {
    medium_computation = ppjsdm::compute_dispersion_for_vcov(medium_dispersion_model, number_types, configuration, nthreads);
  }

  // Finally, add A2 and A3 by using precomputed values above.
  for(std::remove_cv_t<decltype(non_zero_ij_responses.size())> index(0); index < non_zero_ij_responses.size(); ++index) {
    const auto k(non_zero_ij_responses[index]);
    const auto pr(ppjsdm::decode_linear(k, regressors.nrow()));
    const auto i(pr.first);
    const auto j(pr.second);
    const ppjsdm::Marked_point point_i(x[i], y[i], type[i] - 1, mark[i]);
    const ppjsdm::Marked_point point_j(x[j], y[j], type[j] - 1, mark[j]);

    // TODO: Code below is horrible, trying to quickly implement new technique
    size_t i_index_configuration{};
    size_t j_index_configuration{};
    for(size_t l(0); l < ppjsdm::size(configuration); ++l) {
      if(ppjsdm::is_equal(point_i, configuration[l])) {
        i_index_configuration = l;
      }
      if(ppjsdm::is_equal(point_j, configuration[l])) {
        j_index_configuration = l;
      }
    }
    const auto index_in_computation(ppjsdm::encode_linear(std::min(i_index_configuration, j_index_configuration),
                                                          std::max(i_index_configuration, j_index_configuration),
                                                          ppjsdm::size(configuration)));

    using dispersion_t = std::remove_cv_t<std::remove_reference_t<decltype(short_computation.first[0])>>;
    dispersion_t short_i, medium_i, short_j, medium_j;
    if(compute_some_alphas) {
      short_i = short_computation.first[index_in_computation];
      short_j = short_computation.second[index_in_computation];
    }
    if(compute_some_gammas) {
      medium_i = medium_computation.first[index_in_computation];
      medium_j = medium_computation.second[index_in_computation];
    }

    // TODO: Might want to write a function synchronized with prepare_gibbsm_data
    // that constructs the two vectors below.
    std::vector<double> t_i(total_parameters);
    std::vector<double> t_j(total_parameters);

    t_i[ppjsdm::get_type(point_i)] = 1.;
    t_j[ppjsdm::get_type(point_j)] = 1.;

    for(int k2(0); k2 < covariates.size(); ++k2) {
      t_i[index_start_covariates + k2 * number_types + ppjsdm::get_type(point_i)] = covariates_on_configuration[i][k2];
      t_j[index_start_covariates + k2 * number_types + ppjsdm::get_type(point_j)] = covariates_on_configuration[j][k2];
    }

    size_t current_index_alpha(0);
    size_t current_index_gamma(0);
    for(int k1(0); k1 < number_types; ++k1) {
      if(k1 == ppjsdm::get_type(point_i)) {
        size_t filling_alpha(current_index_alpha);
        size_t filling_gamma(current_index_gamma);
        for(int k2(k1); k2 < number_types; ++k2) {
          if(estimate_alpha(k1, k2)) {
            t_i[number_types + filling_alpha] = short_i[k2];
            ++filling_alpha;
          }
          if(estimate_gamma(k1, k2)) {
            t_i[index_start_gamma + filling_gamma] = medium_i[k2];
            ++filling_gamma;
          }
        }
      } else {
        size_t filling_alpha(current_index_alpha);
        size_t filling_gamma(current_index_gamma);
        for(int k2(k1); k2 < number_types; ++k2) {
          if(estimate_alpha(k1, k2)) {
            if(k2 == ppjsdm::get_type(point_i)) {
              t_i[number_types + filling_alpha] = short_i[k1];
            }
            ++filling_alpha;
          }
          if(estimate_gamma(k1, k2)) {
            if(k2 == ppjsdm::get_type(point_i)) {
              t_i[index_start_gamma + filling_gamma] = medium_i[k1];
            }
            ++filling_gamma;
          }
        }
      }

      if(k1 == ppjsdm::get_type(point_j)) {
        for(int k2(k1); k2 < number_types; ++k2) {
          if(estimate_alpha(k1, k2)) {
            t_j[number_types + current_index_alpha] = short_j[k2];
            ++current_index_alpha;
          }
          if(estimate_gamma(k1, k2)) {
            t_j[index_start_gamma + current_index_gamma] = medium_j[k2];
            ++current_index_gamma;
          }
        }
      } else {
        for(int k2(k1); k2 < number_types; ++k2) {
          if(estimate_alpha(k1, k2)) {
            if(k2 == ppjsdm::get_type(point_j)) {
              t_j[number_types + current_index_alpha] = short_j[k1];
            }
            ++current_index_alpha;
          }
          if(estimate_gamma(k1, k2)) {
            if(k2 == ppjsdm::get_type(point_j)) {
              t_j[index_start_gamma + current_index_gamma] = medium_j[k1];
            }
            ++current_index_gamma;
          }
        }
      }
    }

    // TODO: I'm pretty sure papangelou_i_minus_two = exp(t_i, check)...
    double inner_product_i(0);
    double inner_product_j(0);
    for(std::remove_cv_t<decltype(coefficients_vector.size())> k1(0); k1 < coefficients_vector.size(); ++k1) {
      inner_product_i += coefficients_vector[k1] * t_i[k1];
      inner_product_j += coefficients_vector[k1] * t_j[k1];
    }
    const auto papangelou_i(std::exp(inner_product_i));
    const auto papangelou_j(std::exp(inner_product_j));
    const auto constant(2. * (papangelou_i / papangelou[i] - 1.) / ((papangelou_i + rho) * (papangelou_j + rho)));
    for(size_t k1(0); k1 < total_parameters; ++k1) {
      for(size_t k2(0); k2 < total_parameters; ++k2) {
        const auto A2_summand(t_i[k1] * t_j[k2] * constant);
        const auto A3_summand(2. * (t_over_papangelou_on_configuration[i][k1] - t_i[k1] / (papangelou_i + rho)) * (t_over_papangelou_on_configuration[j][k2] - t_j[k2] / (papangelou_j + rho)));
        A(k1, k2) += A2_summand + A3_summand;
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
                                 Rcpp::NumericMatrix short_range,
                                 Rcpp::NumericMatrix medium_range,
                                 Rcpp::NumericMatrix long_range,
                                 R_xlen_t saturation,
                                 double rho,
                                 Rcpp::NumericVector coefficients_vector,
                                 Rcpp::NumericMatrix regressors,
                                 Rcpp::List data_list,
                                 Rcpp::LogicalMatrix estimate_alpha,
                                 Rcpp::LogicalMatrix estimate_gamma,
                                 int nthreads) {
  // Construct std::vector of configurations.
  const ppjsdm::Configuration_wrapper wrapped_configuration(Rcpp::wrap(configuration));
  const auto length_configuration(ppjsdm::size(wrapped_configuration));

  // Convert configuration to std::vector in order for parallelised version to work.
  std::vector<ppjsdm::Marked_point> vector_configuration(length_configuration);
  for(decltype(ppjsdm::size(wrapped_configuration)) j(0); j < length_configuration; ++j) {
    vector_configuration[j] = wrapped_configuration[j];
  }

  const auto number_types(short_range.nrow());
  const auto dispersion(ppjsdm::Saturated_model(model, short_range, saturation));
  const auto medium_range_dispersion(ppjsdm::Saturated_model(medium_range_model, medium_range, long_range, saturation));

  return compute_vcov_helper(vector_configuration,
                             ppjsdm::Im_list_wrapper(covariates),
                             dispersion,
                             medium_range_dispersion,
                             number_types,
                             rho,
                             coefficients_vector,
                             regressors,
                             data_list,
                             estimate_alpha,
                             estimate_gamma,
                             nthreads);
}
