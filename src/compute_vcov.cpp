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

namespace detail {

inline auto make_papangelou(Rcpp::NumericMatrix regressors, Rcpp::NumericVector theta) {
  if(regressors.ncol() != theta.size()) {
    Rcpp::stop("Incompatible sizes when computing the Papangelou intensity");
  }

  std::vector<double> papangelou(regressors.nrow());
  using size_t = decltype(regressors.ncol());
  for(size_t row(0); row < regressors.nrow(); ++row) {
    double inner_product(0);
    for(size_t col(0); col < regressors.ncol(); ++col) {
      inner_product += theta[col] * regressors(row, col);
    }
    papangelou[row] = std::exp(inner_product);
  }

  return papangelou;
}

template<typename Vector>
inline auto make_G2(const Vector& papangelou, Rcpp::NumericVector rho, Rcpp::NumericMatrix regressors, Rcpp::IntegerVector type) {
  const auto number_parameters(regressors.ncol());
  using size_t = decltype(regressors.ncol());

  // G2_temp is a temporary matrix we use to avoid taking the sqrt of rho
  arma::mat G2(number_parameters, number_parameters, arma::fill::zeros);
  arma::mat G2_temp(number_parameters, number_parameters, arma::fill::zeros);

  double kappa(0);
  for(size_t i(0); i < regressors.nrow(); ++i) {
    kappa += 1. / (papangelou[i] + rho[type[i] - 1]);

    const auto constant_1(papangelou[i] / ((papangelou[i] + rho[type[i] - 1]) * (papangelou[i] + rho[type[i] - 1])));
    const auto constant_2(-constant_1 / rho[type[i] - 1]);
    for(size_t k1(0); k1 < number_parameters; ++k1) {
      const auto value_1(regressors(i, k1) * constant_1);
      const auto value_2(regressors(i, k1) * constant_2);
      for(size_t k2(0); k2 < number_parameters; ++k2) {
        G2(k1, k2) += value_2;
        G2_temp(k2, k1) += value_1;
      }
    }
  }

  // The line below should be equivalent to A = -A_temp * A_temp.t() / rho,
  // these complications are to avoid taking the sqrt of rho when it's not uniform.
  G2 *= G2_temp / kappa;

  for(size_t i(0); i < regressors.nrow(); ++i) {
    const auto constant(papangelou[i] * papangelou[i] / ((papangelou[i] + rho[type[i] - 1]) * (papangelou[i] + rho[type[i] - 1]) * (papangelou[i] + rho[type[i] - 1])) / rho[type[i] - 1]);
    for(size_t k1(0); k1 < number_parameters; ++k1) {
      const auto value(regressors(i, k1) * constant);
      for(size_t k2(k1); k2 < number_parameters; ++k2) {
        G2(k1, k2) += value * regressors(i, k2);
      }
    }
  }

  // Symmetrize G2
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    for(size_t k2(k1 + 1); k2 < number_parameters; ++k2) {
      G2(k2, k1) = G2(k1, k2);
    }
  }
  return G2;
}

template<typename Vector>
inline auto make_S(const Vector& papangelou, Rcpp::NumericVector rho, Rcpp::NumericMatrix regressors, Rcpp::IntegerVector type) {
  const auto number_parameters(regressors.ncol());
  using size_t = decltype(regressors.ncol());

  arma::mat S(number_parameters, number_parameters, arma::fill::zeros);

  // Fill upper triangular part of S
  for(size_t row(0); row < regressors.nrow(); ++row) {
    const auto constant(papangelou[row] / ((papangelou[row] + rho[type[row] - 1]) * (papangelou[row] + rho[type[row] - 1])));
    for(size_t k1(0); k1 < number_parameters; ++k1) {
      const auto value(regressors(row, k1) * constant);
      for(size_t k2(k1); k2 < number_parameters; ++k2) {
        S(k1, k2) += value * regressors(row, k2);
      }
    }
  }

  // Symmetrize S
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    for(size_t k2(k1 + 1); k2 < number_parameters; ++k2) {
      S(k2, k1) = S(k1, k2);
    }
  }

  return S;
}

template<typename Vector>
inline auto make_A1(const Vector& papangelou, Rcpp::NumericVector rho, Rcpp::NumericMatrix regressors, Rcpp::IntegerVector type) {
  const auto number_parameters(regressors.ncol());
  using size_t = decltype(regressors.ncol());

  arma::mat A1(number_parameters, number_parameters, arma::fill::zeros);

  // Fill upper triangular part of A1
  for(size_t row(0); row < regressors.nrow(); ++row) {
    const auto constant(papangelou[row] / ((papangelou[row] + rho[type[row] - 1]) * (papangelou[row] + rho[type[row] - 1]) * (papangelou[row] + rho[type[row] - 1])));
    for(size_t k1(0); k1 < number_parameters; ++k1) {
      for(size_t k2(k1); k2 < number_parameters; ++k2) {
        A1(k1, k2) += regressors(row, k1) * regressors(row, k2) * constant;
      }
    }
  }

  // Symmetrize A1
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    for(size_t k2(k1 + 1); k2 < number_parameters; ++k2) {
      A1(k2, k1) = A1(k1, k2);
    }
  }

  return A1;
}

template<typename Vector, typename Configuration>
inline auto make_A2_plus_A3(const Vector& papangelou,
                            Rcpp::NumericVector rho,
                            Rcpp::NumericVector theta,
                            Rcpp::NumericMatrix regressors,
                            Rcpp::List data_list,
                            Rcpp::LogicalMatrix estimate_alpha,
                            Rcpp::LogicalMatrix estimate_gamma,
                            const ppjsdm::Saturated_model& dispersion_model,
                            const ppjsdm::Saturated_model& medium_dispersion_model,
                            const Configuration& configuration,
                            int nthreads,
                            const ppjsdm::Im_list_wrapper& covariates) {
  const auto number_types(rho.size());
  const auto number_parameters_struct(ppjsdm::get_number_parameters(number_types, covariates.size(), estimate_alpha, estimate_gamma));
  const auto index_start_gamma(number_parameters_struct.index_start_gamma);
  const auto index_start_covariates(number_parameters_struct.index_start_covariates);
  const auto number_parameters(number_parameters_struct.total_parameters);

  const auto response(Rcpp::as<Rcpp::IntegerVector>(data_list["response"]));
  const auto x(Rcpp::as<Rcpp::NumericVector>(data_list["x"]));
  const auto y(Rcpp::as<Rcpp::NumericVector>(data_list["y"]));
  const auto type(Rcpp::as<Rcpp::IntegerVector>(data_list["type"]));
  const auto mark(Rcpp::as<Rcpp::NumericVector>(data_list["mark"]));

  using size_t = decltype(ppjsdm::size(configuration));

  arma::mat A(number_parameters, number_parameters, arma::fill::zeros);

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
  using precomputation_t = decltype(ppjsdm::compute_dispersion_for_vcov(dispersion_model, number_types, configuration, nthreads));
  precomputation_t short_computation, medium_computation;
  if(compute_some_alphas) {
    short_computation = ppjsdm::compute_dispersion_for_vcov(dispersion_model, number_types, configuration, nthreads);
  }
  if(compute_some_gammas) {
    medium_computation = ppjsdm::compute_dispersion_for_vcov(medium_dispersion_model, number_types, configuration, nthreads);
  }

  // Finally, add A2 and A3 by using precomputed values above.
  for(decltype(regressors.nrow()) i(0); i < regressors.nrow(); ++i) {
    for(decltype(regressors.nrow()) j(i + 1); j < regressors.nrow(); ++j) {
      if(response[i] == 1 && response[j] == 1) {
        const ppjsdm::Marked_point point_i(x[i], y[i], type[i] - 1, mark[i]);
        const ppjsdm::Marked_point point_j(x[j], y[j], type[j] - 1, mark[j]);

        // Recover the precomputed short and medium range dispersions
        using dispersion_t = std::remove_cv_t<std::remove_reference_t<decltype(short_computation.first[0])>>;
        dispersion_t short_i, medium_i, short_j, medium_j;
        if(compute_some_alphas || compute_some_gammas) {
          size_t i_index_configuration{}, j_index_configuration{};
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
          if(compute_some_alphas) {
            short_i = short_computation.first[index_in_computation];
            short_j = short_computation.second[index_in_computation];
          }
          if(compute_some_gammas) {
            medium_i = medium_computation.first[index_in_computation];
            medium_j = medium_computation.second[index_in_computation];
          }
        }


        // TODO: Might want to write a function synchronized with prepare_gibbsm_data
        // that constructs the two vectors below.
        std::vector<double> t_i(number_parameters);
        std::vector<double> t_j(number_parameters);

        t_i[ppjsdm::get_type(point_i)] = 1.;
        t_j[ppjsdm::get_type(point_j)] = 1.;

        for(int k2(0); k2 < covariates.size(); ++k2) {
          t_i[index_start_covariates + k2 * number_types + ppjsdm::get_type(point_i)] = regressors(i, index_start_covariates + k2 * number_types + ppjsdm::get_type(point_i));
          t_j[index_start_covariates + k2 * number_types + ppjsdm::get_type(point_j)] = regressors(j, index_start_covariates + k2 * number_types + ppjsdm::get_type(point_j));
        }

        size_t current_index_alpha(0);
        size_t current_index_gamma(0);
        for(int k1(0); k1 < number_types; ++k1) {
          if(k1 == ppjsdm::get_type(point_i)) {
            for(int k2(k1); k2 < number_types; ++k2) {
              if(estimate_alpha(k1, k2)) {
                t_i[number_types + current_index_alpha + k2 - k1] = short_i[k2];
              }
              if(estimate_gamma(k1, k2)) {
                t_i[index_start_gamma + current_index_gamma + k2 - k1] = medium_i[k2];
              }
            }
          } else if(k1 < ppjsdm::get_type(point_i)) {
            if(estimate_alpha(k1, ppjsdm::get_type(point_i))) {
              t_i[number_types + current_index_alpha + ppjsdm::get_type(point_i) - k1] = short_i[k1];
            }
            if(estimate_gamma(k1, ppjsdm::get_type(point_i))) {
              t_i[index_start_gamma + current_index_gamma + ppjsdm::get_type(point_i) - k1] = medium_i[k1];
            }
          }

          if(k1 == ppjsdm::get_type(point_j)) {
            for(int k2(k1); k2 < number_types; ++k2) {
              if(estimate_alpha(k1, k2)) {
                t_j[number_types + current_index_alpha + k2 - k1] = short_j[k2];
              }
              if(estimate_gamma(k1, k2)) {
                t_j[index_start_gamma + current_index_gamma + k2 - k1] = medium_j[k2];
              }
            }
          } else if(k1 < ppjsdm::get_type(point_j)) {
            if(estimate_alpha(k1, ppjsdm::get_type(point_j))) {
              t_j[number_types + current_index_alpha + ppjsdm::get_type(point_j) - k1] = short_j[k1];
            }
            if(estimate_gamma(k1, ppjsdm::get_type(point_j))) {
              t_j[index_start_gamma + current_index_gamma + ppjsdm::get_type(point_j) - k1] = medium_j[k1];
            }
          }
          current_index_alpha += number_types - k1;
          current_index_gamma += number_types - k1;
        }

        double inner_product_i(0);
        double inner_product_j(0);
        for(std::remove_cv_t<decltype(theta.size())> k1(0); k1 < theta.size(); ++k1) {
          inner_product_i += theta[k1] * t_i[k1];
          inner_product_j += theta[k1] * t_j[k1];
        }
        const auto papangelou_i(std::exp(inner_product_i));
        const auto papangelou_j(std::exp(inner_product_j));
        const auto constant_1((papangelou_i / papangelou[i] - 1.) / ((papangelou_i + rho[ppjsdm::get_type(point_i)]) * (papangelou_j + rho[ppjsdm::get_type(point_j)])));
        const auto constant_2((papangelou_j / papangelou[j] - 1.) / ((papangelou_j + rho[ppjsdm::get_type(point_j)]) * (papangelou_i + rho[ppjsdm::get_type(point_i)])));
        for(size_t k1(0); k1 < number_parameters; ++k1) {
          for(size_t k2(k1); k2 < number_parameters; ++k2) {
            const auto A2_summand_1(t_i[k1] * t_j[k2] * constant_1);
            const auto A2_summand_2(t_j[k1] * t_i[k2] * constant_2);
            const auto A3_summand_1((regressors(i, k1) / (papangelou[i] + rho[ppjsdm::get_type(point_i)]) - t_i[k1] / (papangelou_i + rho[ppjsdm::get_type(point_i)])) * (regressors(j, k2) / (papangelou[j] + rho[ppjsdm::get_type(point_j)]) - t_j[k2] / (papangelou_j + rho[ppjsdm::get_type(point_j)])));
            const auto A3_summand_2((regressors(j, k1) / (papangelou[j] + rho[ppjsdm::get_type(point_j)]) - t_j[k1] / (papangelou_j + rho[ppjsdm::get_type(point_j)])) * (regressors(i, k2) / (papangelou[i] + rho[ppjsdm::get_type(point_i)]) - t_i[k2] / (papangelou_i + rho[ppjsdm::get_type(point_i)])));
            A(k1, k2) += A2_summand_1 + A2_summand_2 + A3_summand_1 + A3_summand_2;
          }
        }
      }
    }
  }

  // Symmetrize A
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    for(size_t k2(k1 + 1); k2 < number_parameters; ++k2) {
      A(k2, k1) = A(k1, k2);
    }
  }

  return A;
}

} // namespace detail

template<typename Configuration>
Rcpp::NumericMatrix compute_vcov_helper(const Configuration& configuration,
                                        const ppjsdm::Im_list_wrapper& covariates,
                                        const ppjsdm::Saturated_model& dispersion_model,
                                        const ppjsdm::Saturated_model& medium_dispersion_model,
                                        Rcpp::NumericVector rho,
                                        Rcpp::NumericVector theta,
                                        Rcpp::NumericMatrix regressors,
                                        Rcpp::List data_list,
                                        Rcpp::LogicalMatrix estimate_alpha,
                                        Rcpp::LogicalMatrix estimate_gamma,
                                        int nthreads) {
  const auto papangelou(detail::make_papangelou(regressors, theta));
  const auto S(detail::make_S(papangelou, rho, regressors, Rcpp::as<Rcpp::IntegerVector>(data_list["type"])));
  const auto S_inv(arma::inv_sympd(S));
  const auto G2(detail::make_G2(papangelou, rho, regressors, Rcpp::as<Rcpp::IntegerVector>(data_list["type"])));
  const auto A1(detail::make_A1(papangelou, rho, regressors, Rcpp::as<Rcpp::IntegerVector>(data_list["type"])));
  const auto A2_plus_A3(detail::make_A2_plus_A3(papangelou,
                                                rho,
                                                theta,
                                                regressors,
                                                data_list,
                                                estimate_alpha,
                                                estimate_gamma,
                                                dispersion_model,
                                                medium_dispersion_model,
                                                configuration,
                                                nthreads,
                                                covariates));

  Rcpp::NumericMatrix rcpp_matrix(Rcpp::wrap(S_inv * (A1 + A2_plus_A3 + G2) * S_inv));
  const auto col_names(make_model_coloumn_names(covariates, estimate_alpha, estimate_gamma));
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
                                 Rcpp::NumericVector rho,
                                 Rcpp::NumericVector theta,
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

  const auto dispersion(ppjsdm::Saturated_model(model, short_range, saturation));
  const auto medium_range_dispersion(ppjsdm::Saturated_model(medium_range_model, medium_range, long_range, saturation));

  return compute_vcov_helper(vector_configuration,
                             ppjsdm::Im_list_wrapper(covariates),
                             dispersion,
                             medium_range_dispersion,
                             rho,
                             theta,
                             regressors,
                             data_list,
                             estimate_alpha,
                             estimate_gamma,
                             nthreads);
}

