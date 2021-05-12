// [[Rcpp::depends("RcppArmadillo")]]
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

#ifdef _OPENMP
#include <omp.h>
#endif

namespace detail {

template<typename Configuration>
inline void restrict_window(ppjsdm::Window& window, const Configuration& configuration, typename Configuration::size_type max_size, double percent) {
  if(ppjsdm::size(configuration) <= max_size) {
    return;
  }

  while(true) {
    window.shrink_by_percent(percent);
    decltype(max_size) restricted_size(0);
    for(const auto& point: configuration) {
      if(window.is_in(ppjsdm::get_x(point), ppjsdm::get_y(point))) {
        ++restricted_size;
        if(restricted_size > max_size) {
          break;
        }
      }
    }
    if(restricted_size <= max_size) {
      return;
    }
  }
}

template<typename Configuration>
inline auto restrict_configuration(const ppjsdm::Window& restricted_window, const Configuration& configuration) {
  // Create object to return and reserve
  Configuration restricted_configuration;
  restricted_configuration.reserve(ppjsdm::size(configuration));

  // Fill in with points within the window
  for(const auto& point: configuration) {
    if(restricted_window.is_in(ppjsdm::get_x(point), ppjsdm::get_y(point))) {
      restricted_configuration.emplace_back(point);
    }
  }

  return restricted_configuration;
}

inline auto make_papangelou(const ppjsdm::Lightweight_matrix<double>& regressors,
                            const std::vector<double>& theta,
                            int nthreads) {
  if(regressors.ncol() != static_cast<decltype(regressors.ncol())>(theta.size())) {
    Rcpp::stop("Incompatible sizes when computing the Papangelou intensity");
  }

#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif

  using computation_t = long double;
  using size_t = decltype(regressors.ncol());

  std::vector<double> papangelou(regressors.nrow());

#pragma omp parallel
{
  decltype(papangelou) papangelou_private(regressors.nrow());
#pragma omp for nowait
  for(size_t row = 0; row < regressors.nrow(); ++row) {
    computation_t inner_product(0.);
    for(size_t col(0); col < regressors.ncol(); ++col) {
      const auto a = static_cast<computation_t>(theta[col]);
      const auto b = static_cast<computation_t>(regressors(row, col));
      inner_product += a * b;
    }
    papangelou_private[row] = static_cast<typename decltype(papangelou)::value_type>(std::exp(inner_product));
  }
#pragma omp critical
  for(size_t row(0); row < regressors.nrow(); ++row) {
    papangelou[row] += papangelou_private[row];
  }
}

  return papangelou;
}

inline auto make_G2(const std::vector<double>& papangelou,
                    const std::vector<double>& rho,
                    const ppjsdm::Lightweight_matrix<double>& regressors,
                    const std::vector<int>& type,
                    double window_volume,
                    int nthreads) {
  using computation_t = long double;
  using size_t = decltype(regressors.ncol());

#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif

  const auto number_parameters(regressors.ncol());

  std::vector<computation_t> G2_vec(number_parameters);
  computation_t kappa(0.);

  // Note: Strangely enough the two parallel sections below cannot be combined.
  // I'm getting unexplained bugs involving NaN values when I try to.
#pragma omp parallel
{
  decltype(G2_vec) G2_vec_private(number_parameters);
  decltype(kappa) kappa_private(0.);
#pragma omp for nowait
  for(size_t i = 0; i < regressors.nrow(); ++i) {
    const auto papangelou_value = static_cast<computation_t>(papangelou[i]);
    const auto papangelou_plus_rho = papangelou_value + static_cast<computation_t>(rho[type[i] - 1]);
    const auto ratio(papangelou_value / papangelou_plus_rho);
    kappa_private += static_cast<computation_t>(rho[type[i] - 1]) * ratio;

    const auto constant = ratio / papangelou_plus_rho;
    for(size_t k1(0); k1 < number_parameters; ++k1) {
      G2_vec_private[k1] += static_cast<computation_t>(regressors(i, k1)) * constant;
    }
  }
#pragma omp critical
  kappa += kappa_private;
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    G2_vec[k1] += G2_vec_private[k1];
  }
}

  // TODO: Square matrix instead
  ppjsdm::Lightweight_matrix<computation_t> G2(number_parameters, number_parameters);

#pragma omp parallel
{
  decltype(G2) G2_private(number_parameters, number_parameters);
#pragma omp for nowait
  for(size_t i = 0; i < regressors.nrow(); ++i) {
    const auto papangelou_value = static_cast<computation_t>(papangelou[i]);
    const auto papangelou_plus_rho = papangelou_value + static_cast<computation_t>(rho[type[i] - 1]);
    const auto ratio(papangelou_value / papangelou_plus_rho);
    const auto constant(ratio * ratio / papangelou_plus_rho / static_cast<computation_t>(rho[type[i] - 1]));
    for(size_t k1(0); k1 < number_parameters; ++k1) {
      const auto value = static_cast<computation_t>(regressors(i, k1)) * constant;
      for(size_t k2(k1); k2 < number_parameters; ++k2) {
        G2_private(k1, k2) += value * static_cast<computation_t>(regressors(i, k2));
      }
    }
  }
#pragma omp critical
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    for(size_t k2(k1); k2 < number_parameters; ++k2) {
      G2(k1, k2) += G2_private(k1, k2);
    }
  }
}

  const auto a = static_cast<computation_t>(window_volume) * static_cast<computation_t>(number_parameters) / kappa;
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    const auto b = a * G2_vec[k1];
    for(size_t k2(k1); k2 < number_parameters; ++k2) {
      G2(k1, k2) -= b * G2_vec[k2];
    }
  }

  // Construct the return object
  arma::mat A(number_parameters, number_parameters);
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    A(k1, k1) = static_cast<typename decltype(A)::value_type>(G2(k1, k1));
    for(size_t k2(k1 + 1); k2 < number_parameters; ++k2) {
      A(k1, k2) = static_cast<typename decltype(A)::value_type>(G2(k1, k2));
      A(k2, k1) = A(k1, k2);
    }
  }

  return A;
}

inline auto make_S(const std::vector<double>& papangelou,
                   const std::vector<double>& rho,
                   const ppjsdm::Lightweight_matrix<double>& regressors,
                   const std::vector<int>& type,
                   int nthreads) {
  using computation_t = long double;
  using size_t = decltype(regressors.ncol());

#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif

  const auto number_parameters(regressors.ncol());

  ppjsdm::Lightweight_matrix<computation_t> S(number_parameters, number_parameters);

#pragma omp parallel
{
  decltype(S) S_private(number_parameters, number_parameters);
#pragma omp for nowait
  for(size_t row = 0; row < regressors.nrow(); ++row) {
    const auto papangelou_value(static_cast<computation_t>(papangelou[row]));
    const auto papangelou_plus_rho(papangelou_value + static_cast<computation_t>(rho[type[row] - 1]));
    const auto constant(papangelou_value / (papangelou_plus_rho * papangelou_plus_rho));
    for(size_t k1(0); k1 < number_parameters; ++k1) {
      const auto value(static_cast<computation_t>(regressors(row, k1)) * constant);
      for(size_t k2(k1); k2 < number_parameters; ++k2) {
        S_private(k1, k2) += value * static_cast<computation_t>(regressors(row, k2));
      }
    }
  }
#pragma omp critical
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    for(size_t k2(k1); k2 < number_parameters; ++k2) {
      S(k1, k2) += S_private(k1, k2);
    }
  }
}

  // Construct the return object
  arma::mat A(number_parameters, number_parameters);
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    A(k1, k1) = static_cast<typename decltype(A)::value_type>(S(k1, k1));
    for(size_t k2(k1 + 1); k2 < number_parameters; ++k2) {
      A(k1, k2) = static_cast<typename decltype(A)::value_type>(S(k1, k2));
      A(k2, k1) = A(k1, k2);
    }
  }

  return A;
}

inline auto make_A1(const std::vector<double>& papangelou,
                    const std::vector<double>& rho,
                    const ppjsdm::Lightweight_matrix<double>& regressors,
                    const std::vector<int>& type,
                    int nthreads) {
  using computation_t = long double;
  using size_t = decltype(regressors.ncol());

#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif

  const auto number_parameters(regressors.ncol());

  ppjsdm::Lightweight_matrix<computation_t> A1(number_parameters, number_parameters);

#pragma omp parallel
{
  decltype(A1) A1_private(number_parameters, number_parameters);
#pragma omp for nowait
  for(size_t row = 0; row < regressors.nrow(); ++row) {
    const auto papangelou_value = static_cast<computation_t>(papangelou[row]);
    const auto papangelou_plus_rho = papangelou_value + static_cast<computation_t>(rho[type[row] - 1]);
    const auto constant(papangelou_value / (papangelou_plus_rho * papangelou_plus_rho * papangelou_plus_rho));
    for(size_t k1(0); k1 < number_parameters; ++k1) {
      const auto value = static_cast<computation_t>(regressors(row, k1) * constant);
      for(size_t k2(k1); k2 < number_parameters; ++k2) {
        A1_private(k1, k2) += value * static_cast<computation_t>(regressors(row, k2));
      }
    }
  }
#pragma omp critical
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    for(size_t k2(k1); k2 < number_parameters; ++k2) {
      A1(k1, k2) += A1_private(k1, k2);
    }
  }
}

  // Construct the return object
  arma::mat A(number_parameters, number_parameters);
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    A(k1, k1) = static_cast<typename decltype(A)::value_type>(A1(k1, k1));
    for(size_t k2(k1 + 1); k2 < number_parameters; ++k2) {
      A(k1, k2) = static_cast<typename decltype(A)::value_type>(A1(k1, k2));
      A(k2, k1) = A(k1, k2);
    }
  }

  return A;
}

template<typename Configuration>
inline auto make_A2_plus_A3(const std::vector<double>& papangelou,
                            const std::vector<double>& rho,
                            const std::vector<double>& theta,
                            const ppjsdm::Lightweight_matrix<double>& regressors,
                            Rcpp::List data_list,
                            const ppjsdm::Lightweight_matrix<bool>& estimate_alpha,
                            const ppjsdm::Lightweight_matrix<bool>& estimate_gamma,
                            const ppjsdm::Saturated_model& dispersion_model,
                            const ppjsdm::Saturated_model& medium_dispersion_model,
                            const Configuration& configuration,
                            int nthreads,
                            const ppjsdm::Im_list_wrapper& covariates,
                            double initial_window_volume,
                            const ppjsdm::Window& restricted_window) {
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif

  const auto restricted_configuration(detail::restrict_configuration(restricted_window, configuration));

  const int number_types(rho.size());
  const auto number_parameters_struct(ppjsdm::get_number_parameters(number_types, covariates.size(), estimate_alpha, estimate_gamma));
  const auto index_start_gamma(number_parameters_struct.index_start_gamma);
  const auto index_start_covariates(number_parameters_struct.index_start_covariates);
  const auto number_parameters(number_parameters_struct.total_parameters);

  const auto response(Rcpp::as<std::vector<int>>(data_list["response"]));
  const auto x(Rcpp::as<std::vector<double>>(data_list["x"]));
  const auto y(Rcpp::as<std::vector<double>>(data_list["y"]));
  const auto type(Rcpp::as<std::vector<int>>(data_list["type"]));
  const auto mark(Rcpp::as<std::vector<double>>(data_list["mark"]));

  using computation_t = long double;
  using size_t = decltype(ppjsdm::size(configuration));

  // Automatically initialized to zeros
  ppjsdm::Lightweight_matrix<computation_t> mat(number_parameters, number_parameters);

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

  // Add A2 and A3 by using precomputed values above.
  // TODO: This depends on memory
  // TODO: Clean up code below, new technique
  const long long int increment(1000000);
  for(long long int filling(0); filling < static_cast<long long int>(ppjsdm::size(restricted_configuration)) * (static_cast<long long int>(ppjsdm::size(restricted_configuration)) - 1) / 2; filling += increment) {
    // Parallel computations done before adding A2 and A3
    using precomputation_t = decltype(ppjsdm::compute_dispersion_for_vcov(dispersion_model, number_types, configuration, nthreads));
    precomputation_t short_computation, medium_computation;
    if(compute_some_alphas) {
      short_computation = ppjsdm::compute_dispersion_for_vcov(dispersion_model, number_types, configuration, restricted_configuration, filling, filling + increment, nthreads);
    }
    if(compute_some_gammas) {
      medium_computation = ppjsdm::compute_dispersion_for_vcov(medium_dispersion_model, number_types, configuration, restricted_configuration, filling, filling + increment, nthreads);
    }

    const auto max_index(std::min<long long int>(filling + increment, static_cast<long long int>(ppjsdm::size(restricted_configuration)) * (static_cast<long long int>(ppjsdm::size(restricted_configuration)) - 1) / 2));
#pragma omp parallel
{
    decltype(mat) mat_private(number_parameters, number_parameters);
#pragma omp for nowait
    for(long long int index_in_computation = filling; index_in_computation < max_index; ++index_in_computation) {
      const auto pr(ppjsdm::decode_linear(index_in_computation, ppjsdm::size(restricted_configuration)));
      const size_t i(pr.first);
      const size_t j(pr.second);
      const auto point_i(restricted_configuration[i]);
      const auto point_j(restricted_configuration[j]);
      long long int i_index_regressors{};
      for(decltype(regressors.nrow()) l(0); l < regressors.nrow(); ++l) {
        if(ppjsdm::is_equal(point_i, ppjsdm::Marked_point(x[l], y[l], type[l] - 1, mark[l]))) {
          i_index_regressors = l;
          break;
        }
      }
      long long int j_index_regressors{};
      for(decltype(regressors.nrow()) l(0); l < regressors.nrow(); ++l) {
        if(ppjsdm::is_equal(point_j, ppjsdm::Marked_point(x[l], y[l], type[l] - 1, mark[l]))) {
          j_index_regressors = l;
          break;
        }
      }

      // Recover the precomputed short and medium range dispersions
      using dispersion_t = std::remove_cv_t<std::remove_reference_t<decltype(short_computation.first[0])>>;
      dispersion_t short_i, medium_i, short_j, medium_j;
      if(compute_some_alphas) {
        short_i = short_computation.first[index_in_computation - filling];
        short_j = short_computation.second[index_in_computation - filling];
      }
      if(compute_some_gammas) {
        medium_i = medium_computation.first[index_in_computation - filling];
        medium_j = medium_computation.second[index_in_computation - filling];
      }


      // TODO: Might want to write a function synchronized with prepare_gibbsm_data
      // that constructs the two vectors below.
      std::vector<computation_t> t_i(number_parameters);
      std::vector<computation_t> t_j(number_parameters);

      t_i[ppjsdm::get_type(point_i)] = static_cast<computation_t>(1.);
      t_j[ppjsdm::get_type(point_j)] = static_cast<computation_t>(1.);

      for(int k2(0); k2 < covariates.size(); ++k2) {
        t_i[index_start_covariates + k2 * number_types + ppjsdm::get_type(point_i)] = static_cast<computation_t>(covariates[k2](point_i));
        t_j[index_start_covariates + k2 * number_types + ppjsdm::get_type(point_j)] = static_cast<computation_t>(covariates[k2](point_j));
      }

      size_t current_index_alpha(0);
      size_t current_index_gamma(0);
      for(int k1(0); k1 < number_types; ++k1) {
        if(k1 == ppjsdm::get_type(point_i)) {
          for(int k2(k1); k2 < number_types; ++k2) {
            if(estimate_alpha(k1, k2)) {
              t_i[number_types + current_index_alpha + k2 - k1] = static_cast<computation_t>(short_i[k2]);
            }
            if(estimate_gamma(k1, k2)) {
              t_i[index_start_gamma + current_index_gamma + k2 - k1] = static_cast<computation_t>(medium_i[k2]);
            }
          }
        } else if(k1 < ppjsdm::get_type(point_i)) {
          if(estimate_alpha(k1, ppjsdm::get_type(point_i))) {
            t_i[number_types + current_index_alpha + ppjsdm::get_type(point_i) - k1] = static_cast<computation_t>(short_i[k1]);
          }
          if(estimate_gamma(k1, ppjsdm::get_type(point_i))) {
            t_i[index_start_gamma + current_index_gamma + ppjsdm::get_type(point_i) - k1] = static_cast<computation_t>(medium_i[k1]);
          }
        }

        if(k1 == ppjsdm::get_type(point_j)) {
          for(int k2(k1); k2 < number_types; ++k2) {
            if(estimate_alpha(k1, k2)) {
              t_j[number_types + current_index_alpha + k2 - k1] = static_cast<computation_t>(short_j[k2]);
            }
            if(estimate_gamma(k1, k2)) {
              t_j[index_start_gamma + current_index_gamma + k2 - k1] = static_cast<computation_t>(medium_j[k2]);
            }
          }
        } else if(k1 < ppjsdm::get_type(point_j)) {
          if(estimate_alpha(k1, ppjsdm::get_type(point_j))) {
            t_j[number_types + current_index_alpha + ppjsdm::get_type(point_j) - k1] = static_cast<computation_t>(short_j[k1]);
          }
          if(estimate_gamma(k1, ppjsdm::get_type(point_j))) {
            t_j[index_start_gamma + current_index_gamma + ppjsdm::get_type(point_j) - k1] = static_cast<computation_t>(medium_j[k1]);
          }
        }
        current_index_alpha += number_types - k1;
        current_index_gamma += number_types - k1;
      }

      computation_t inner_product_i(0);
      computation_t inner_product_j(0);
      for(std::remove_cv_t<decltype(theta.size())> k1(0); k1 < theta.size(); ++k1) {
        inner_product_i += static_cast<computation_t>(theta[k1]) * t_i[k1];
        inner_product_j += static_cast<computation_t>(theta[k1]) * t_j[k1];
      }
      const auto papangelou_i(std::exp(inner_product_i));
      const auto papangelou_j(std::exp(inner_product_j));
      const auto a1 = papangelou_i / static_cast<computation_t>(papangelou[i_index_regressors]) - static_cast<computation_t>(1.);
      const auto a2 = (papangelou_i + static_cast<computation_t>(rho[ppjsdm::get_type(point_i)])) * (papangelou_j + static_cast<computation_t>(rho[ppjsdm::get_type(point_j)]));
      const auto a3 = papangelou_j / static_cast<computation_t>(papangelou[j_index_regressors]) - static_cast<computation_t>(1.);
      const auto constant_1 = a1 / a2;
      const auto constant_2 = a3 / a2;
      for(size_t k1(0); k1 < number_parameters; ++k1) {
        const auto a1 = static_cast<computation_t>(regressors(i_index_regressors, k1)) / (static_cast<computation_t>(papangelou[i_index_regressors]) + static_cast<computation_t>(rho[ppjsdm::get_type(point_i)]));
        const auto a2 = t_i[k1] / (papangelou_i + static_cast<computation_t>(rho[ppjsdm::get_type(point_i)]));
        const auto c1 = static_cast<computation_t>(regressors(j_index_regressors, k1)) / (static_cast<computation_t>(papangelou[j_index_regressors]) + static_cast<computation_t>(rho[ppjsdm::get_type(point_j)]));
        const auto c2 = t_j[k1] / (papangelou_j + static_cast<computation_t>(rho[ppjsdm::get_type(point_j)]));
        for(size_t k2(k1); k2 < number_parameters; ++k2) {
          mat_private(k1, k2) += t_i[k1] * t_j[k2] * constant_1;
          mat_private(k1, k2) += t_j[k1] * t_i[k2] * constant_2;
          const auto b1 = static_cast<computation_t>(regressors(j_index_regressors, k2)) / (static_cast<computation_t>(papangelou[j_index_regressors]) + static_cast<computation_t>(rho[ppjsdm::get_type(point_j)]));
          const auto b2 = t_j[k2] / (papangelou_j + static_cast<computation_t>(rho[ppjsdm::get_type(point_j)]));
          mat_private(k1, k2) += (a1 - a2) * (b1 - b2);
          const auto d1 = static_cast<computation_t>(regressors(i_index_regressors, k2)) / (static_cast<computation_t>(papangelou[i_index_regressors]) + static_cast<computation_t>(rho[ppjsdm::get_type(point_i)]));
          const auto d2 = t_i[k2] / (papangelou_i + static_cast<computation_t>(rho[ppjsdm::get_type(point_i)]));
          mat_private(k1, k2) += (c1 - c2) * (d1 - d2);
        }
      }
    }
#pragma omp critical
    for(size_t k1(0); k1 < number_parameters; ++k1) {
      for(size_t k2(k1); k2 < number_parameters; ++k2) {
        mat(k1, k2) += mat_private(k1, k2);
      }
    }
}
  }

  arma::mat A(number_parameters, number_parameters);

  // Fill and symmetrize A
  const auto normalisation_constant(initial_window_volume / restricted_window.volume());
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    A(k1, k1) = normalisation_constant * mat(k1, k1);
    for(size_t k2(k1 + 1); k2 < number_parameters; ++k2) {
      A(k1, k2) = normalisation_constant * mat(k1, k2);
      A(k2, k1) = A(k1, k2);
    }
  }

  return A;
}

} // namespace detail

template<typename Configuration>
Rcpp::List compute_vcov_helper(const Configuration& configuration,
                               ppjsdm::Window& window,
                               const ppjsdm::Im_list_wrapper& covariates,
                               const ppjsdm::Saturated_model& dispersion_model,
                               const ppjsdm::Saturated_model& medium_dispersion_model,
                               const std::vector<double>& rho,
                               double initial_window_volume,
                               const std::vector<double>& theta,
                               const ppjsdm::Lightweight_matrix<double>& regressors,
                               Rcpp::List data_list,
                               const ppjsdm::Lightweight_matrix<bool>& estimate_alpha,
                               const ppjsdm::Lightweight_matrix<bool>& estimate_gamma,
                               int nthreads) {
  const auto papangelou(detail::make_papangelou(regressors, theta, nthreads));
  const auto S(detail::make_S(papangelou, rho, regressors, Rcpp::as<std::vector<int>>(data_list["type"]), nthreads));
  if(S.has_nan()) {
    Rcpp::stop("Found NaN values in matrix S (in vcov computation).");
  }
  const auto S_inv(arma::inv_sympd(S));
  const auto G2(detail::make_G2(papangelou, rho, regressors, Rcpp::as<std::vector<int>>(data_list["type"]), initial_window_volume, nthreads));
  const auto A1(detail::make_A1(papangelou, rho, regressors, Rcpp::as<std::vector<int>>(data_list["type"]), nthreads));

  detail::restrict_window(window, configuration, 1000, 0.05);
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
                                                covariates,
                                                initial_window_volume,
                                                window));

  const auto col_names(make_model_coloumn_names(covariates, estimate_alpha, estimate_gamma));

  Rcpp::NumericMatrix G1_rcpp(Rcpp::wrap(S_inv * (A1 + A2_plus_A3) * S_inv));
  Rcpp::NumericMatrix G2_rcpp(Rcpp::wrap(S_inv * G2 * S_inv));

  Rcpp::colnames(G1_rcpp) = col_names;
  Rcpp::rownames(G1_rcpp) = col_names;
  Rcpp::colnames(G2_rcpp) = col_names;
  Rcpp::rownames(G2_rcpp) = col_names;

  return Rcpp::List::create(Rcpp::Named("G1") = G1_rcpp,
                            Rcpp::Named("G2") = G2_rcpp);
}

// [[Rcpp::export]]
Rcpp::List compute_vcov(SEXP configuration,
                        SEXP window,
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
  // Convert the SEXP configuration to a C++ object.
  const ppjsdm::Configuration_wrapper wrapped_configuration(Rcpp::wrap(configuration));
  const auto length_configuration(ppjsdm::size(wrapped_configuration));

  // Convert configuration to std::vector in order for parallelised version to work.
  std::vector<ppjsdm::Marked_point> vector_configuration(length_configuration);
  for(decltype(ppjsdm::size(wrapped_configuration)) j(0); j < length_configuration; ++j) {
    vector_configuration[j] = wrapped_configuration[j];
  }

  // Convert the SEXP window to a C++ object.
  ppjsdm::Window cpp_window(window);

  // Call the main function.
  return compute_vcov_helper(vector_configuration,
                             cpp_window,
                             ppjsdm::Im_list_wrapper(covariates),
                             ppjsdm::Saturated_model(model, short_range, saturation),
                             ppjsdm::Saturated_model(medium_range_model, medium_range, long_range, saturation),
                             Rcpp::as<std::vector<double>>(rho),
                             cpp_window.volume(),
                             Rcpp::as<std::vector<double>>(theta),
                             ppjsdm::Lightweight_matrix<double>(regressors),
                             data_list,
                             ppjsdm::Lightweight_matrix<bool>(estimate_alpha),
                             ppjsdm::Lightweight_matrix<bool>(estimate_gamma),
                             nthreads);
}

