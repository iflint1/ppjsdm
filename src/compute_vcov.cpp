// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rinternals.h>

#include "configuration/configuration_manipulation.hpp"
#include "configuration/configuration_wrapper.hpp"
#include "configuration/get_number_points.hpp"

#include "point/point_manipulation.hpp"

#include "saturated_model/compute_dispersion_vcov.hpp"
#include "saturated_model/exponential_family_model.hpp"
#include "saturated_model/regression_vector.hpp"
#include "saturated_model/saturated_model.hpp"

#include "simulation/rstrat_single.hpp"

#include "utility/construct_if_missing.hpp"
#include "utility/flatten_strict_upper_triangular.hpp"
#include "utility/gibbsm_helper_functions.hpp"
#include "utility/im_wrapper.hpp"
#include "utility/is_symmetric_matrix.hpp"
#include "utility/lightweight_matrix.hpp"
#include "utility/subset_to_nonNA.hpp"
#include "utility/timer.hpp"
#include "utility/window.hpp"
#include "utility/window_concrete.hpp"

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

inline auto make_dispersions(Rcpp::List model,
                             Rcpp::List short_range,
                             R_xlen_t saturation) {
  using dispersion_t = ppjsdm::Saturated_model<double>;
  std::vector<dispersion_t> dispersions;
  for(decltype(short_range.size()) i(0); i < short_range.size(); ++i) {
    dispersions.emplace_back(dispersion_t(Rcpp::as<Rcpp::CharacterVector>(model[i]),
                                          Rcpp::as<Rcpp::NumericMatrix>(short_range[i]),
                                          saturation));
  }
  return dispersions;
}

inline auto convert_list_of_boolean_matrices_to_vector(Rcpp::List boolean_matrices) {
  using conversion_t = ppjsdm::Lightweight_matrix<bool>;
  std::vector<conversion_t> converted_matrices(boolean_matrices.size());
  for(decltype(converted_matrices.size()) i(0); i < converted_matrices.size(); ++i) {
    converted_matrices[i] = conversion_t(Rcpp::as<Rcpp::LogicalMatrix>(boolean_matrices[i]));
  }
  return converted_matrices;
}

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

template<typename Configuration, typename Window>
inline auto restrict_configuration(const Window& restricted_window, const Configuration& configuration) {
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

  using computation_t = double;
  using size_t = decltype(regressors.ncol());

  std::vector<double> papangelou(regressors.nrow());

#pragma omp parallel default(none) shared(regressors, theta, papangelou)
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

inline auto make_G2_binomial(const std::vector<double>& papangelou,
                             const std::vector<double>& rho,
                             const ppjsdm::Lightweight_matrix<double>& regressors,
                             const std::vector<int>& type,
                             double window_volume,
                             int nthreads) {
  using computation_t = double;
  using size_t = decltype(regressors.ncol());

#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif

  const auto number_parameters(regressors.ncol());

  std::vector<computation_t> G2_vec(number_parameters);
  ppjsdm::Lightweight_square_matrix<computation_t> G2(number_parameters);
  computation_t kappa(0.);

#pragma omp parallel shared(papangelou, rho, regressors, type, kappa, G2, G2_vec)
{
  decltype(G2_vec) G2_vec_private(number_parameters);
  decltype(G2) G2_private(number_parameters);
  decltype(kappa) kappa_private(0.);
#pragma omp for nowait
  for(size_t i = 0; i < regressors.nrow(); ++i) {
    const auto papangelou_value = static_cast<computation_t>(papangelou[i]);
    const auto papangelou_plus_rho = papangelou_value + static_cast<computation_t>(rho[type[i] - 1]);
    const auto ratio(papangelou_value / papangelou_plus_rho);
    kappa_private += static_cast<computation_t>(rho[type[i] - 1]) * ratio;

    const auto constant(ratio * ratio / papangelou_plus_rho / static_cast<computation_t>(rho[type[i] - 1]));
    const auto constant_vec = ratio / papangelou_plus_rho;
    for(size_t k1(0); k1 < number_parameters; ++k1) {
      G2_vec_private[k1] += static_cast<computation_t>(regressors(i, k1)) * constant_vec;

      const auto value = static_cast<computation_t>(regressors(i, k1)) * constant;
      for(size_t k2(k1); k2 < number_parameters; ++k2) {
        G2_private(k1, k2) += value * static_cast<computation_t>(regressors(i, k2));
      }
    }
  }
#pragma omp critical
{
  kappa += kappa_private;
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    G2_vec[k1] += G2_vec_private[k1];
    for(size_t k2(k1); k2 < number_parameters; ++k2) {
      G2(k1, k2) += G2_private(k1, k2);
    }
  }
}
}

const auto a = static_cast<computation_t>(number_parameters) / kappa / static_cast<computation_t>(window_volume);
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

template<typename Configuration, typename FloatType, typename VectorBoolMatrices>
inline auto make_G2_stratified(const Configuration& configuration,
                               const Configuration& dummy,
                               const ppjsdm::Window& window,
                               const std::vector<double>& theta,
                               const std::vector<double>& rho,
                               const ppjsdm::Lightweight_matrix<double>& regressors,
                               const VectorBoolMatrices& estimate_alpha,
                               const ppjsdm::Lightweight_square_matrix<bool>& estimate_gamma,
                               const std::vector<ppjsdm::Saturated_model<FloatType>>& dispersion_model,
                               const ppjsdm::Saturated_model<FloatType>& medium_dispersion_model,
                               const ppjsdm::Im_list_wrapper& covariates,
                               int nthreads) {
  using computation_t = double;
  using size_t = decltype(regressors.ncol());

  // Extract some values relating to the number of parameters and how they're ordered
  const int number_types(rho.size());
  const auto number_parameters_struct(ppjsdm::get_number_parameters(number_types, covariates.size(), estimate_alpha, estimate_gamma));
  const auto index_start_gamma(number_parameters_struct.index_start_gamma);
  const auto index_start_covariates(number_parameters_struct.index_start_covariates);
  const auto number_parameters(number_parameters_struct.total_parameters);

#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif

  // Construct other stratified point process
  // std::vector<computation_t> delta(rho.size());
  // for(decltype(rho.size()) i(0); i < rho.size(); ++i) {
  //   delta[i] = static_cast<computation_t>(1.) / std::sqrt(static_cast<computation_t>(rho[i]));
  // }
  const auto other_stratified(subset_to_nonNA(ppjsdm::rstratpp_single<Configuration>(window, ppjsdm::get_number_points(dummy)), covariates));

  if(ppjsdm::size(dummy) != ppjsdm::size(other_stratified)) {
    Rcpp::Rcout << "Size dummy: " << ppjsdm::size(dummy) << " and size new stratified: " << ppjsdm::size(other_stratified) << ".\n";
    Rcpp::stop("The dummy points and the independent draw of a stratified binomial point process should have the same number of points. This could be due to NA values on the covariates, causing some of the dummy points to be removed on the draws. To solve this, either subset the window to locations where the covariates are non-NA, or avoid dummy points distributed as a stratified point process.");
  }

  ppjsdm::Lightweight_square_matrix<computation_t> G2(number_parameters);

  // Do you need to compute some of the alphas/gammas? If not, we can skip some of the dispersion computations
  bool compute_some_alphas(false);
  bool compute_some_gammas(false);
  for(int i(0); i < number_types; ++i) {
    for(int j(i); j < number_types; ++j) {
      for(decltype(estimate_alpha.size()) k(0); k < estimate_alpha.size(); ++k) {
        if(estimate_alpha[k](i, j)) {
          compute_some_alphas = true;
        }
      }
      if(estimate_gamma(i, j)) {
        compute_some_gammas = true;
      }
    }
  }

  // Precompute dispersions
  using precomputation_t = decltype(ppjsdm::compute_dispersion_for_fitting<false>(medium_dispersion_model, 1, 1, configuration, dummy));
  std::vector<precomputation_t> short_computation_dummy, short_computation_other;
  precomputation_t medium_computation_dummy, medium_computation_other;

  if(compute_some_alphas) {
    // TODO: compute_some_alphas can be made to be specific to k
    for(decltype(estimate_alpha.size()) k(0); k < estimate_alpha.size(); ++k) {
      short_computation_dummy.emplace_back(ppjsdm::compute_dispersion_for_fitting<false>(dispersion_model[k], number_types, nthreads, configuration, dummy));
      short_computation_other.emplace_back(ppjsdm::compute_dispersion_for_fitting<false>(dispersion_model[k], number_types, nthreads, configuration, other_stratified));
    }
  }
  if(compute_some_gammas) {
    medium_computation_dummy = ppjsdm::compute_dispersion_for_fitting<false>(medium_dispersion_model, number_types, nthreads, configuration, dummy);
    medium_computation_other = ppjsdm::compute_dispersion_for_fitting<false>(medium_dispersion_model, number_types, nthreads, configuration, other_stratified);
  }

  // Compute G2
#pragma omp parallel shared(G2, dummy, short_computation_dummy, short_computation_other) \
  shared(medium_computation_dummy, medium_computation_other, compute_some_alphas, compute_some_gammas) \
  shared(covariates, estimate_alpha, estimate_gamma, theta, rho)
{
  decltype(G2) G2_private(number_parameters);
#pragma omp for nowait
  for(decltype(ppjsdm::size(dummy)) i = 0; i < ppjsdm::size(dummy); ++i) {
    // Recover the precomputed short and medium range dispersions
    using dispersion_t = std::remove_cv_t<std::remove_reference_t<decltype(medium_computation_dummy[0])>>;
    std::vector<dispersion_t> short_dummy, short_other;
    dispersion_t medium_dummy, medium_other;
    if(compute_some_alphas) {
      for(decltype(estimate_alpha.size()) k(0); k < estimate_alpha.size(); ++k) {
        short_dummy.emplace_back(short_computation_dummy[k][i]);
        short_other.emplace_back(short_computation_other[k][i]);
      }
    }
    if(compute_some_gammas) {
      medium_dummy = medium_computation_dummy[i];
      medium_other = medium_computation_other[i];
    }

    std::vector<double> dummy_covariates(covariates.size());
    std::vector<double> other_covariates(covariates.size());
    for(decltype(covariates.size()) k2(0); k2 < covariates.size(); ++k2) {
      dummy_covariates[k2] = covariates[k2](dummy[i]);
      other_covariates[k2] = covariates[k2](other_stratified[i]);
    }

    const auto t_dummy(ppjsdm::make_regression_vector<computation_t>(ppjsdm::get_type(dummy[i]),
                                                                     number_types,
                                                                     number_parameters,
                                                                     index_start_covariates,
                                                                     index_start_gamma,
                                                                     dummy_covariates,
                                                                     short_dummy,
                                                                     medium_dummy,
                                                                     estimate_alpha,
                                                                     estimate_gamma));
    const auto t_other(ppjsdm::make_regression_vector<computation_t>(ppjsdm::get_type(other_stratified[i]),
                                                                     number_types,
                                                                     number_parameters,
                                                                     index_start_covariates,
                                                                     index_start_gamma,
                                                                     other_covariates,
                                                                     short_other,
                                                                     medium_other,
                                                                     estimate_alpha,
                                                                     estimate_gamma));

    computation_t inner_product_dummy(0);
    computation_t inner_product_other(0);
    for(size_t k1(0); k1 < number_parameters; ++k1) {
      inner_product_dummy += static_cast<computation_t>(theta[k1]) * t_dummy[k1];
      inner_product_other += static_cast<computation_t>(theta[k1]) * t_other[k1];
    }
    const auto papangelou_dummy(std::exp(inner_product_dummy));
    const auto papangelou_other(std::exp(inner_product_other));
    const auto papangelou_dummy_plus_rho(papangelou_dummy + static_cast<computation_t>(rho[ppjsdm::get_type(dummy[i])]));
    const auto papangelou_other_plus_rho(papangelou_other + static_cast<computation_t>(rho[ppjsdm::get_type(other_stratified[i])]));

    std::vector<computation_t> values(number_parameters);
    for(size_t k1(0); k1 < number_parameters; ++k1) {
      values[k1] = t_dummy[k1] * papangelou_dummy / papangelou_dummy_plus_rho - t_other[k1] * papangelou_other / papangelou_other_plus_rho;
    }

    const auto one_over_rho_squared(static_cast<computation_t>(1.) / static_cast<computation_t>(rho[ppjsdm::get_type(dummy[i])]) / static_cast<computation_t>(rho[ppjsdm::get_type(dummy[i])]));
    for(size_t k1(0); k1 < number_parameters; ++k1) {
      for(size_t k2(k1); k2 < number_parameters; ++k2) {
        G2_private(k1, k2) += values[k1] * values[k2] * one_over_rho_squared;
      }
    }
  }
#pragma omp critical
{
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    for(size_t k2(k1); k2 < number_parameters; ++k2) {
      G2(k1, k2) += G2_private(k1, k2);
    }
  }
}
}

  // Construct the return object
  arma::mat A(number_parameters, number_parameters);
  for(size_t k1(0); k1 < number_parameters; ++k1) {
    A(k1, k1) = static_cast<typename decltype(A)::value_type>(static_cast<computation_t>(0.5) * G2(k1, k1));
    for(size_t k2(k1 + 1); k2 < number_parameters; ++k2) {
      A(k1, k2) = static_cast<typename decltype(A)::value_type>(static_cast<computation_t>(0.5) * G2(k1, k2));
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
  using computation_t = double;
  using size_t = decltype(regressors.ncol());

#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif

  const auto number_parameters(regressors.ncol());

  ppjsdm::Lightweight_square_matrix<computation_t> S(number_parameters);

#pragma omp parallel shared(papangelou, rho, regressors, type, S)
{
  decltype(S) S_private(number_parameters);
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
  using computation_t = double;
  using size_t = decltype(regressors.ncol());

#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif

  const auto number_parameters(regressors.ncol());

  ppjsdm::Lightweight_square_matrix<computation_t> A1(number_parameters);

#pragma omp parallel shared(papangelou, rho, regressors, type, A1)
{
  decltype(A1) A1_private(number_parameters);
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

template<typename Configuration, typename FloatType, typename Window, typename VectorBoolMatrices>
inline auto make_A2_plus_A3(const std::vector<double>& papangelou,
                            const std::vector<double>& rho,
                            const std::vector<double>& theta,
                            const ppjsdm::Lightweight_matrix<double>& regressors,
                            Rcpp::List data_list,
                            const VectorBoolMatrices& estimate_alpha,
                            const ppjsdm::Lightweight_square_matrix<bool>& estimate_gamma,
                            const std::vector<ppjsdm::Saturated_model<FloatType>>& dispersion_model,
                            const ppjsdm::Saturated_model<FloatType>& medium_dispersion_model,
                            const Configuration& configuration,
                            int nthreads,
                            const ppjsdm::Im_list_wrapper& covariates,
                            double initial_window_volume,
                            const Window& restricted_window,
                            bool debug) {
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif

  // Timer used to print debugging information
  ppjsdm::PreciseTimer timer{};

  // A2 and A3 are going to be computed with only points from the restricted window.
  const auto restricted_configuration(detail::restrict_configuration(restricted_window, configuration));

  // Extract some values relating to the number of parameters and how they're ordered
  const int number_types(rho.size());
  const auto number_parameters_struct(ppjsdm::get_number_parameters(number_types, covariates.size(), estimate_alpha, estimate_gamma));
  const auto index_start_gamma(number_parameters_struct.index_start_gamma);
  const auto index_start_covariates(number_parameters_struct.index_start_covariates);
  const auto number_parameters(number_parameters_struct.total_parameters);

  // We're given an R List containing the following vectors, convert them to thread-safe C++ containers
  const auto response(Rcpp::as<std::vector<int>>(data_list["response"]));
  const auto x(Rcpp::as<std::vector<double>>(data_list["x"]));
  const auto y(Rcpp::as<std::vector<double>>(data_list["y"]));
  const auto type(Rcpp::as<std::vector<int>>(data_list["type"]));
  const auto mark(Rcpp::as<std::vector<double>>(data_list["mark"]));

  // A few useful typedefs
  using computation_t = double;
  using size_t = decltype(ppjsdm::size(configuration));

  // Do we need to compute some of the alphas/gammas? If not, we can skip some of the dispersion computations
  bool compute_some_alphas(false);
  bool compute_some_gammas(false);
  for(int i(0); i < number_types; ++i) {
    for(int j(i); j < number_types; ++j) {
      for(decltype(estimate_alpha.size()) k(0); k < estimate_alpha.size(); ++k) {
        if(estimate_alpha[k](i, j)) {
          compute_some_alphas = true;
        }
      }
      if(estimate_gamma(i, j)) {
        compute_some_gammas = true;
      }
    }
  }

  // In the main loop, we need to compute papangelou[i] and regressors[i, *] / (rho + papangelou[i])
  // for indices i that correspond to points in the restricted configuration. Precompute this here.
  std::vector<std::vector<computation_t>> regressors_over_papangelou_plus_rho(ppjsdm::size(restricted_configuration));
  std::vector<computation_t> restricted_papangelou(ppjsdm::size(restricted_configuration));
  for(size_t i(0); i < ppjsdm::size(restricted_configuration); ++i) {
    regressors_over_papangelou_plus_rho[i] = typename decltype(regressors_over_papangelou_plus_rho)::value_type(number_parameters);

    using x_size_t = typename decltype(x)::size_type;
    x_size_t index(-1);
    for(x_size_t l(0); l < x.size(); ++l) {
      if(ppjsdm::is_equal(restricted_configuration[i], ppjsdm::Marked_point(x[l], y[l], type[l] - 1, mark[l]))) {
        index = l;
        break;
      }
    }
    if(index == static_cast<x_size_t>(-1)) {
      Rcpp::Rcout << index << '\n';
      Rcpp::stop("Did not find an element of the restricted configuration in the complete configuration.");
    }

    const auto r = static_cast<computation_t>(rho[ppjsdm::get_type(restricted_configuration[i])]);
    restricted_papangelou[i] = static_cast<computation_t>(papangelou[index]);
    for(size_t k1(0); k1 < number_parameters; ++k1) {
      const auto reg = static_cast<computation_t>(regressors(index, k1));
      regressors_over_papangelou_plus_rho[i][k1] = reg / (restricted_papangelou[i] + r);
    }
  }

  // Automatically initialized to zeros
  ppjsdm::Lightweight_square_matrix<computation_t> mat(number_parameters);

  // Add A2 and A3 by using precomputed values above.
  // TODO: This depends on memory
  // TODO: Clean up code below, new technique
  const long long int increment(10000000);
  for(long long int filling(0); filling < static_cast<long long int>(ppjsdm::size(restricted_configuration)) * (static_cast<long long int>(ppjsdm::size(restricted_configuration)) - 1) / 2; filling += increment) {
    // Parallel computations done before adding A2 and A3
    using precomputation_t = decltype(ppjsdm::compute_dispersion_for_vcov(medium_dispersion_model, number_types, configuration, restricted_configuration, 0, 0, 1));
    std::vector<precomputation_t> short_computation;
    precomputation_t medium_computation;

    if(debug) {
      timer.set_current();
      Rcpp::Rcout << "Starting computation of batch of dispersions... Number of points in batch: " << ppjsdm::size(restricted_configuration) << ".\n";
    }
    if(compute_some_alphas) {
      for(decltype(estimate_alpha.size()) k(0); k < estimate_alpha.size(); ++k) {
        short_computation.emplace_back(ppjsdm::compute_dispersion_for_vcov(dispersion_model[k], number_types, configuration, restricted_configuration, filling, filling + increment, nthreads, debug));
      }
    }
    if(compute_some_gammas) {
      medium_computation = ppjsdm::compute_dispersion_for_vcov(medium_dispersion_model, number_types, configuration, restricted_configuration, filling, filling + increment, nthreads, debug);
    }
    if(debug) {
      Rcpp::Rcout << "Computed batch of dispersions. Elapsed time: " << timer.print_elapsed_time();
      timer.set_current();
      Rcpp::Rcout << "Filling the matrix...\n";
    }

    const auto max_index(std::min<long long int>(filling + increment, static_cast<long long int>(ppjsdm::size(restricted_configuration)) * (static_cast<long long int>(ppjsdm::size(restricted_configuration)) - 1) / 2));
#pragma omp parallel                                              \
    shared(filling, regressors, compute_some_alphas, compute_some_gammas)           \
      shared(short_computation, medium_computation, estimate_alpha, estimate_gamma) \
      shared(covariates, theta, papangelou, rho, mat, restricted_papangelou, regressors_over_papangelou_plus_rho)
{
        decltype(mat) mat_private(number_parameters);
#pragma omp for nowait
        for(long long int index_in_computation = filling; index_in_computation < max_index; ++index_in_computation) {
          const auto pr(ppjsdm::decode_linear(index_in_computation, ppjsdm::size(restricted_configuration)));
          const size_t i(pr.first);
          const size_t j(pr.second);
          const auto point_i(restricted_configuration[i]);
          const auto point_j(restricted_configuration[j]);

          // Recover the precomputed short and medium range dispersions
          using dispersion_t = std::remove_cv_t<std::remove_reference_t<decltype(medium_computation.first[0])>>;
          std::vector<dispersion_t> short_i, short_j;
          dispersion_t medium_i, medium_j;
          if(compute_some_alphas) {
            for(decltype(estimate_alpha.size()) k(0); k < estimate_alpha.size(); ++k) {
              short_i.emplace_back(short_computation[k].first[index_in_computation - filling]);
              short_j.emplace_back(short_computation[k].second[index_in_computation - filling]);
            }
          }
          if(compute_some_gammas) {
            medium_i = medium_computation.first[index_in_computation - filling];
            medium_j = medium_computation.second[index_in_computation - filling];
          }

          std::vector<double> covariates_i(covariates.size());
          std::vector<double> covariates_j(covariates.size());
          for(decltype(covariates.size()) k2(0); k2 < covariates.size(); ++k2) {
            covariates_i[k2] = covariates[k2](point_i);
            covariates_j[k2] = covariates[k2](point_j);
          }

          const auto t_i(ppjsdm::make_regression_vector<computation_t>(ppjsdm::get_type(point_i),
                                                                       number_types,
                                                                       number_parameters,
                                                                       index_start_covariates,
                                                                       index_start_gamma,
                                                                       covariates_i,
                                                                       short_i,
                                                                       medium_i,
                                                                       estimate_alpha,
                                                                       estimate_gamma));
          const auto t_j(ppjsdm::make_regression_vector<computation_t>(ppjsdm::get_type(point_j),
                                                                       number_types,
                                                                       number_parameters,
                                                                       index_start_covariates,
                                                                       index_start_gamma,
                                                                       covariates_j,
                                                                       short_j,
                                                                       medium_j,
                                                                       estimate_alpha,
                                                                       estimate_gamma));

          computation_t inner_product_i(0);
          computation_t inner_product_j(0);
          for(size_t k1(0); k1 < number_parameters; ++k1) {
            inner_product_i += static_cast<computation_t>(theta[k1]) * t_i[k1];
            inner_product_j += static_cast<computation_t>(theta[k1]) * t_j[k1];
          }
          const auto papangelou_i(std::exp(inner_product_i));
          const auto papangelou_j(std::exp(inner_product_j));
          const auto papangelou_i_plus_rho(papangelou_i + static_cast<computation_t>(rho[ppjsdm::get_type(point_i)]));
          const auto papangelou_j_plus_rho(papangelou_j + static_cast<computation_t>(rho[ppjsdm::get_type(point_j)]));

          // Precompute delta
          std::vector<computation_t> delta_i(number_parameters), delta_j(number_parameters);
          for(size_t k1(0); k1 < number_parameters; ++k1) {
            delta_i[k1] = regressors_over_papangelou_plus_rho[i][k1] - t_i[k1] / papangelou_i_plus_rho;
            delta_j[k1] = regressors_over_papangelou_plus_rho[j][k1] - t_j[k1] / papangelou_j_plus_rho;
          }

          const auto x1 = papangelou_i / restricted_papangelou[i] - static_cast<computation_t>(1.);
          const auto x2 = papangelou_i_plus_rho * papangelou_j_plus_rho;
          const auto constant = x1 / x2;
          for(size_t k1(0); k1 < number_parameters; ++k1) {
            for(size_t k2(k1); k2 < number_parameters; ++k2) {
              mat_private(k1, k2) += (t_i[k1] * t_j[k2] + t_j[k1] * t_i[k2]) * constant;
              mat_private(k1, k2) += delta_i[k1] * delta_j[k2] + delta_j[k1] * delta_i[k2];
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
    if(debug) {
      Rcpp::Rcout << "Finished filling matrix. Elapsed time: " << timer.print_elapsed_time();
      timer.set_current();
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

template<typename Configuration, typename FloatType, typename VectorBoolMatrices>
Rcpp::List compute_vcov_helper(const Configuration& configuration,
                               const Configuration& dummy,
                               ppjsdm::Window& window,
                               const ppjsdm::Im_list_wrapper& covariates,
                               const std::vector<ppjsdm::Saturated_model<FloatType>>& dispersion_model,
                               const ppjsdm::Saturated_model<FloatType>& medium_dispersion_model,
                               const std::vector<double>& rho,
                               double initial_window_volume,
                               const std::vector<double>& theta,
                               const ppjsdm::Lightweight_matrix<double>& regressors,
                               Rcpp::List data_list,
                               const VectorBoolMatrices& estimate_alpha,
                               const ppjsdm::Lightweight_square_matrix<bool>& estimate_gamma,
                               bool debug,
                               int nthreads,
                               int npoints,
                               bool multiple_windows,
                               std::string dummy_distribution,
                               Rcpp::NumericVector mark_range) {
  // Extract some values relating to the number of parameters and how they're ordered
  const auto number_parameters(ppjsdm::get_number_parameters(rho.size(), covariates.size(), estimate_alpha, estimate_gamma).total_parameters);

  ppjsdm::PreciseTimer timer{};
  if(debug) {
    Rcpp::Rcout << "Starting computation of the Papangelou conditional intensity...\n";
  }
  const auto papangelou(detail::make_papangelou(regressors, theta, nthreads));
  if(debug) {
    Rcpp::Rcout << "Finished computing the Papangelou conditional intensity. Time elapsed: " << timer.print_elapsed_time();
    timer.set_current();
    Rcpp::Rcout << "Starting computation of S...\n";
  }
  const auto S(detail::make_S(papangelou, rho, regressors, Rcpp::as<std::vector<int>>(data_list["type"]), nthreads));
  if(debug) {
    Rcpp::Rcout << "Finished computing S. Time elapsed: " << timer.print_elapsed_time();
    timer.set_current();
    Rcpp::Rcout << "Starting inversion of S...\n";
  }
  if(S.has_nan()) {
    Rcpp::stop("Found NaN values in matrix S (in vcov computation).");
  }
  // arma::inv appears to be less risky than arma::inv_sympd,
  // especially in corner cases where S appears numerically to not be positive definite.
  const auto S_inv(arma::inv(S));
  if(debug) {
    Rcpp::Rcout << "Finished inversion of S. Time elapsed: " << timer.print_elapsed_time();
    timer.set_current();
    Rcpp::Rcout << "Starting computation of G2...\n";
  }
  using g2_t = decltype(detail::make_G2_binomial(papangelou, rho, regressors, std::vector<int>{}, 1., 1));
  g2_t G2;
  if(dummy_distribution == std::string("binomial")) {
    G2 = detail::make_G2_binomial(papangelou, rho, regressors, Rcpp::as<std::vector<int>>(data_list["type"]), initial_window_volume, nthreads);
  } else if(dummy_distribution == std::string("stratified")) {
    G2 = detail::make_G2_stratified(configuration, dummy, window, theta, rho, regressors, estimate_alpha, estimate_gamma, dispersion_model, medium_dispersion_model, covariates, nthreads);
  } else {
    Rcpp::stop("Unknown dummy_distribution parameter.");
  }
  if(debug) {
    Rcpp::Rcout << "Finished computation of G2. Time elapsed: " << timer.print_elapsed_time();
    timer.set_current();
    Rcpp::Rcout << "Starting computation of A1...\n";
  }
  const auto A1(detail::make_A1(papangelou, rho, regressors, Rcpp::as<std::vector<int>>(data_list["type"]), nthreads));
  if(debug) {
    Rcpp::Rcout << "Finished computation of A1. Time elapsed: " << timer.print_elapsed_time();
    timer.set_current();
    Rcpp::Rcout << "Starting computation of A2 + A3...\n";
  }

  arma::mat A2_plus_A3;
  if(multiple_windows) {
    A2_plus_A3 = arma::mat(number_parameters, number_parameters, arma::fill::zeros);

    const auto delta_x(window.xmax() - window.xmin());
    const auto delta_y(window.ymax() - window.ymin());
    int nx, ny;
    if(delta_x >= delta_y) {
      const auto r(delta_y / delta_x);
      nx = std::max<int>(1, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(ppjsdm::size(configuration)) / (r * static_cast<double>(npoints))))));
      ny = static_cast<int>(std::ceil(r * nx));
    } else {
      const auto r(delta_x / delta_y);
      ny = std::max<int>(1, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(ppjsdm::size(configuration)) / (r * static_cast<double>(npoints))))));
      nx = static_cast<int>(std::ceil(r * ny));
    }

    for(int i(0); i < nx; ++i) {
      for(int j(0); j < ny; ++j) {
        const ppjsdm::detail::Rectangle_window restricted_window(window.xmin() + static_cast<double>(i) * delta_x / static_cast<double>(nx),
                                                                 window.xmin() + static_cast<double>(i + 1) * delta_x / static_cast<double>(nx),
                                                                 window.ymin() + static_cast<double>(j) * delta_y / static_cast<double>(ny),
                                                                 window.ymin() + static_cast<double>(j + 1) * delta_y / static_cast<double>(ny),
                                                                 mark_range);
        A2_plus_A3 += detail::make_A2_plus_A3(papangelou,
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
                                              restricted_window,
                                              debug);
      }
    }
    A2_plus_A3 /= static_cast<double>(nx * ny);
  } else {
    detail::restrict_window(window, configuration, npoints, 0.05);
    A2_plus_A3 = detail::make_A2_plus_A3(papangelou,
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
                                         window,
                                         debug);
  }

  if(debug) {
    Rcpp::Rcout << "Finished computation of A2 + A3. Time elapsed: " << timer.print_elapsed_time();
    timer.set_current();
    Rcpp::Rcout << "Starting clean-up...\n";
  }

  const auto col_names(make_model_coloumn_names(covariates, estimate_alpha, estimate_gamma));

  Rcpp::NumericMatrix G1_rcpp(Rcpp::wrap(S_inv * (A1 + A2_plus_A3) * S_inv));
  Rcpp::NumericMatrix G2_rcpp(Rcpp::wrap(S_inv * G2 * S_inv));

  Rcpp::colnames(G1_rcpp) = col_names;
  Rcpp::rownames(G1_rcpp) = col_names;
  Rcpp::colnames(G2_rcpp) = col_names;
  Rcpp::rownames(G2_rcpp) = col_names;

  if(debug) {
    Rcpp::Rcout << "Finished computing the vcov matrix. Time elapsed: " << timer.print_elapsed_time();
    Rcpp::Rcout << "Total time taken to compute vcov matrix: " << timer.print_total_time();
  }

  return Rcpp::List::create(Rcpp::Named("G1") = G1_rcpp,
                            Rcpp::Named("G2") = G2_rcpp);
}

// template<typename Configuration, typename FloatType, typename VectorBoolMatrices>
// Rcpp::List compute_S_helper(const Configuration& configuration,
//                             const Configuration& dummy,
//                             ppjsdm::Window& window,
//                             const ppjsdm::Im_list_wrapper& covariates,
//                             const std::vector<ppjsdm::Saturated_model<FloatType>>& dispersion_model,
//                             const ppjsdm::Saturated_model<FloatType>& medium_dispersion_model,
//                             const std::vector<double>& rho,
//                             double initial_window_volume,
//                             const std::vector<double>& theta,
//                             const ppjsdm::Lightweight_matrix<double>& regressors,
//                             Rcpp::List data_list,
//                             const VectorBoolMatrices& estimate_alpha,
//                             const ppjsdm::Lightweight_square_matrix<bool>& estimate_gamma,
//                             bool debug,
//                             int nthreads,
//                             int npoints,
//                             bool multiple_windows,
//                             std::string dummy_distribution,
//                             Rcpp::NumericVector mark_range) {
//   // Extract some values relating to the number of parameters and how they're ordered
//   const auto number_parameters(ppjsdm::get_number_parameters(rho.size(), covariates.size(), estimate_alpha, estimate_gamma).total_parameters);
//
//   ppjsdm::PreciseTimer timer{};
//   if(debug) {
//     Rcpp::Rcout << "Starting computation of the Papangelou conditional intensity...\n";
//   }
//   const auto papangelou(detail::make_papangelou(regressors, theta, nthreads));
//   if(debug) {
//     Rcpp::Rcout << "Finished computing the Papangelou conditional intensity. Time elapsed: " << timer.print_elapsed_time();
//     timer.set_current();
//     Rcpp::Rcout << "Starting computation of S...\n";
//   }
//   const auto S(detail::make_S(papangelou, rho, regressors, Rcpp::as<std::vector<int>>(data_list["type"]), nthreads));
//   if(debug) {
//     Rcpp::Rcout << "Finished computing S. Time elapsed: " << timer.print_elapsed_time();
//     timer.set_current();
//     Rcpp::Rcout << "Starting inversion of S...\n";
//   }
//   if(S.has_nan()) {
//     Rcpp::stop("Found NaN values in matrix S (in vcov computation).");
//   }
//   // arma::inv appears to be less risky than arma::inv_sympd,
//   // especially in corner cases where S appears numerically to not be positive definite.
//   const auto S_inv(arma::inv(S));
//   if(debug) {
//     Rcpp::Rcout << "Finished inversion of S. Time elapsed: " << timer.print_elapsed_time();
//     timer.set_current();
//     Rcpp::Rcout << "Starting computation of G2...\n";
//   }
//   using g2_t = decltype(detail::make_G2_binomial(papangelou, rho, regressors, std::vector<int>{}, 1., 1));
//   g2_t G2;
//   if(dummy_distribution == std::string("binomial")) {
//     G2 = detail::make_G2_binomial(papangelou, rho, regressors, Rcpp::as<std::vector<int>>(data_list["type"]), initial_window_volume, nthreads);
//   } else if(dummy_distribution == std::string("stratified")) {
//     G2 = detail::make_G2_stratified(configuration, dummy, window, theta, rho, regressors, estimate_alpha, estimate_gamma, dispersion_model, medium_dispersion_model, covariates, nthreads);
//   } else {
//     Rcpp::stop("Unknown dummy_distribution parameter.");
//   }
//   if(debug) {
//     Rcpp::Rcout << "Finished computation of G2. Time elapsed: " << timer.print_elapsed_time();
//     timer.set_current();
//     Rcpp::Rcout << "Starting computation of A1...\n";
//   }
//   const auto A1(detail::make_A1(papangelou, rho, regressors, Rcpp::as<std::vector<int>>(data_list["type"]), nthreads));
//   if(debug) {
//     Rcpp::Rcout << "Finished computation of A1. Time elapsed: " << timer.print_elapsed_time();
//     timer.set_current();
//     Rcpp::Rcout << "Starting computation of A2 + A3...\n";
//   }
//
//   arma::mat A2_plus_A3;
//   if(multiple_windows) {
//     A2_plus_A3 = arma::mat(number_parameters, number_parameters, arma::fill::zeros);
//
//     const auto delta_x(window.xmax() - window.xmin());
//     const auto delta_y(window.ymax() - window.ymin());
//     int nx, ny;
//     if(delta_x >= delta_y) {
//       const auto r(delta_y / delta_x);
//       nx = std::max<int>(1, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(ppjsdm::size(configuration)) / (r * static_cast<double>(npoints))))));
//       ny = static_cast<int>(std::ceil(r * nx));
//     } else {
//       const auto r(delta_x / delta_y);
//       ny = std::max<int>(1, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(ppjsdm::size(configuration)) / (r * static_cast<double>(npoints))))));
//       nx = static_cast<int>(std::ceil(r * ny));
//     }
//
//     for(int i(0); i < nx; ++i) {
//       for(int j(0); j < ny; ++j) {
//         const ppjsdm::detail::Rectangle_window restricted_window(window.xmin() + static_cast<double>(i) * delta_x / static_cast<double>(nx),
//                                                                  window.xmin() + static_cast<double>(i + 1) * delta_x / static_cast<double>(nx),
//                                                                  window.ymin() + static_cast<double>(j) * delta_y / static_cast<double>(ny),
//                                                                  window.ymin() + static_cast<double>(j + 1) * delta_y / static_cast<double>(ny),
//                                                                  mark_range);
//         A2_plus_A3 += detail::make_A2_plus_A3(papangelou,
//                                               rho,
//                                               theta,
//                                               regressors,
//                                               data_list,
//                                               estimate_alpha,
//                                               estimate_gamma,
//                                               dispersion_model,
//                                               medium_dispersion_model,
//                                               configuration,
//                                               nthreads,
//                                               covariates,
//                                               initial_window_volume,
//                                               restricted_window,
//                                               debug);
//       }
//     }
//     A2_plus_A3 /= (nx * ny);
//   } else {
//     detail::restrict_window(window, configuration, npoints, 0.05);
//     A2_plus_A3 = detail::make_A2_plus_A3(papangelou,
//                                          rho,
//                                          theta,
//                                          regressors,
//                                          data_list,
//                                          estimate_alpha,
//                                          estimate_gamma,
//                                          dispersion_model,
//                                          medium_dispersion_model,
//                                          configuration,
//                                          nthreads,
//                                          covariates,
//                                          initial_window_volume,
//                                          window,
//                                          debug);
//   }
//
//   if(debug) {
//     Rcpp::Rcout << "Finished computation of A2 + A3. Time elapsed: " << timer.print_elapsed_time();
//     timer.set_current();
//     Rcpp::Rcout << "Starting clean-up...\n";
//   }
//
//   const auto col_names(make_model_coloumn_names(covariates, estimate_alpha, estimate_gamma));
//
//   Rcpp::NumericMatrix G1_rcpp(Rcpp::wrap(S_inv * (A1 + A2_plus_A3) * S_inv));
//   Rcpp::NumericMatrix G2_rcpp(Rcpp::wrap(S_inv * G2 * S_inv));
//
//   Rcpp::colnames(G1_rcpp) = col_names;
//   Rcpp::rownames(G1_rcpp) = col_names;
//   Rcpp::colnames(G2_rcpp) = col_names;
//   Rcpp::rownames(G2_rcpp) = col_names;
//
//   if(debug) {
//     Rcpp::Rcout << "Finished computing the vcov matrix. Time elapsed: " << timer.print_elapsed_time();
//     Rcpp::Rcout << "Total time taken to compute vcov matrix: " << timer.print_total_time();
//   }
//
//   return Rcpp::List::create(Rcpp::Named("G1") = G1_rcpp,
//                             Rcpp::Named("G2") = G2_rcpp);
// }

// [[Rcpp::export]]
Rcpp::List compute_vcov(SEXP configuration,
                        SEXP dummy,
                        SEXP window,
                        Rcpp::List covariates,
                        SEXP model,
                        Rcpp::CharacterVector medium_range_model,
                        Rcpp::List short_range,
                        Rcpp::NumericMatrix medium_range,
                        Rcpp::NumericMatrix long_range,
                        R_xlen_t saturation,
                        Rcpp::NumericVector rho,
                        Rcpp::NumericVector theta,
                        Rcpp::NumericMatrix regressors,
                        Rcpp::List data_list,
                        Rcpp::List estimate_alpha,
                        Rcpp::LogicalMatrix estimate_gamma,
                        bool debug,
                        int nthreads,
                        int npoints,
                        bool multiple_windows,
                        std::string dummy_distribution,
                        Rcpp::NumericVector mark_range) {
  // Convert the SEXP configuration to a C++ object.
  const ppjsdm::Configuration_wrapper wrapped_configuration(Rcpp::wrap(configuration));

  // Convert configuration to std::vector in order for parallelised version to work.
  const auto length_configuration(ppjsdm::size(wrapped_configuration));
  std::vector<ppjsdm::Marked_point> vector_configuration(length_configuration);
  for(decltype(ppjsdm::size(wrapped_configuration)) j(0); j < length_configuration; ++j) {
    vector_configuration[j] = wrapped_configuration[j];
  }

  // Convert the SEXP dummy to a C++ object.
  const ppjsdm::Configuration_wrapper wrapped_dummy(Rcpp::wrap(dummy));

  // Convert configuration to std::vector in order for parallelised version to work.
  const auto length_dummy(ppjsdm::size(wrapped_dummy));
  std::vector<ppjsdm::Marked_point> vector_dummy(length_dummy);
  for(decltype(ppjsdm::size(wrapped_dummy)) j(0); j < length_dummy; ++j) {
    vector_dummy[j] = wrapped_dummy[j];
  }

  // Convert the SEXP window to a C++ object.
  ppjsdm::Window cpp_window(window, mark_range);

  // Call the main function.
  return compute_vcov_helper(vector_configuration,
                             vector_dummy,
                             cpp_window,
                             ppjsdm::Im_list_wrapper(covariates),
                             detail::make_dispersions(model, short_range, saturation),
                             ppjsdm::Saturated_model<double>(medium_range_model, medium_range, long_range, saturation),
                             Rcpp::as<std::vector<double>>(rho),
                             cpp_window.volume(),
                             Rcpp::as<std::vector<double>>(theta),
                             ppjsdm::Lightweight_matrix<double>(regressors),
                             data_list,
                             detail::convert_list_of_boolean_matrices_to_vector(estimate_alpha),
                             ppjsdm::Lightweight_square_matrix<bool>(estimate_gamma),
                             debug,
                             nthreads,
                             npoints,
                             multiple_windows,
                             dummy_distribution,
                             mark_range);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_S_cpp(Rcpp::NumericVector rho,
                                  Rcpp::NumericVector theta,
                                  Rcpp::NumericMatrix regressors,
                                  Rcpp::IntegerVector type,
                                  int nthreads) {
  // Construct papangelou intensity
  const auto papangelou(detail::make_papangelou(ppjsdm::Lightweight_matrix<double>(regressors),
                                                Rcpp::as<std::vector<double>>(theta), nthreads));

  // Compute and return S
  const auto S(detail::make_S(papangelou,
                              Rcpp::as<std::vector<double>>(rho),
                              ppjsdm::Lightweight_matrix<double>(regressors),
                              Rcpp::as<std::vector<int>>(type),
                              nthreads));

  return Rcpp::wrap(S);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_A1_cpp(Rcpp::NumericVector rho,
                                   Rcpp::NumericVector theta,
                                   Rcpp::NumericMatrix regressors,
                                   Rcpp::IntegerVector type,
                                   int nthreads) {
  // Construct papangelou intensity
  const auto papangelou(detail::make_papangelou(ppjsdm::Lightweight_matrix<double>(regressors),
                                                Rcpp::as<std::vector<double>>(theta), nthreads));

  // Compute and return A1
  const auto A1(detail::make_A1(papangelou,
                                Rcpp::as<std::vector<double>>(rho),
                                ppjsdm::Lightweight_matrix<double>(regressors),
                                Rcpp::as<std::vector<int>>(type),
                                nthreads));

  return Rcpp::wrap(A1);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_A2_plus_A3_cpp(SEXP configuration,
                                           SEXP window,
                                           Rcpp::List covariates,
                                           SEXP model,
                                           Rcpp::CharacterVector medium_range_model,
                                           Rcpp::List short_range,
                                           Rcpp::NumericMatrix medium_range,
                                           Rcpp::NumericMatrix long_range,
                                           R_xlen_t saturation,
                                           Rcpp::NumericVector rho,
                                           Rcpp::NumericVector theta,
                                           Rcpp::NumericMatrix regressors,
                                           Rcpp::List data_list,
                                           Rcpp::List estimate_alpha,
                                           Rcpp::LogicalMatrix estimate_gamma,
                                           int nthreads,
                                           int npoints,
                                           bool multiple_windows,
                                           Rcpp::NumericVector mark_range,
                                           bool debug,
                                           int max_executions) {
  const auto number_parameters(ppjsdm::get_number_parameters(rho.size(), covariates.size(),
                                                             detail::convert_list_of_boolean_matrices_to_vector(estimate_alpha),
                                                             estimate_gamma).total_parameters);

  // Convert the SEXP configuration to a C++ object.
  const ppjsdm::Configuration_wrapper wrapped_configuration(Rcpp::wrap(configuration));

  // Convert configuration to std::vector in order for parallelised version to work.
  const auto length_configuration(ppjsdm::size(wrapped_configuration));
  std::vector<ppjsdm::Marked_point> vector_configuration(length_configuration);
  for(decltype(ppjsdm::size(wrapped_configuration)) j(0); j < length_configuration; ++j) {
    vector_configuration[j] = wrapped_configuration[j];
  }

  // Convert the SEXP window to a C++ object.
  ppjsdm::Window cpp_window(window, mark_range);

  const auto initial_window_volume(cpp_window.volume());

  // Construct papangelou intensity
  const auto papangelou(detail::make_papangelou(ppjsdm::Lightweight_matrix<double>(regressors), Rcpp::as<std::vector<double>>(theta), nthreads));

  // Compute and return A2 + A3
  arma::mat A2_plus_A3;
  if(multiple_windows) {
    A2_plus_A3 = arma::mat(number_parameters, number_parameters, arma::fill::zeros);

    const auto delta_x(cpp_window.xmax() - cpp_window.xmin());
    const auto delta_y(cpp_window.ymax() - cpp_window.ymin());
    int nx, ny;
    if(delta_x >= delta_y) {
      const auto r(delta_y / delta_x);
      nx = std::max<int>(1, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(ppjsdm::size(vector_configuration)) / (r * static_cast<double>(npoints))))));
      ny = static_cast<int>(std::ceil(r * nx));
    } else {
      const auto r(delta_x / delta_y);
      ny = std::max<int>(1, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(ppjsdm::size(vector_configuration)) / (r * static_cast<double>(npoints))))));
      nx = static_cast<int>(std::ceil(r * ny));
    }
    int total_executions(0);
    const auto random_rows(Rcpp::sample(nx, nx, false, R_NilValue, false));
    const auto random_cols(Rcpp::sample(ny, ny, false, R_NilValue, false));
    for(int i(0); i < nx; ++i) {
      for(int j(0); j < ny; ++j) {
        if(total_executions < max_executions) {
          const auto random_row(random_rows[i]);
          const auto random_col(random_cols[j]);
          const ppjsdm::detail::Rectangle_window restricted_window(cpp_window.xmin() + static_cast<double>(random_row) * delta_x / static_cast<double>(nx),
                                                                   cpp_window.xmin() + static_cast<double>(random_row + 1) * delta_x / static_cast<double>(nx),
                                                                   cpp_window.ymin() + static_cast<double>(random_col) * delta_y / static_cast<double>(ny),
                                                                   cpp_window.ymin() + static_cast<double>(random_col + 1) * delta_y / static_cast<double>(ny),
                                                                   mark_range);
          A2_plus_A3 += detail::make_A2_plus_A3(papangelou,
                                                Rcpp::as<std::vector<double>>(rho),
                                                Rcpp::as<std::vector<double>>(theta),
                                                ppjsdm::Lightweight_matrix<double>(regressors),
                                                data_list,
                                                detail::convert_list_of_boolean_matrices_to_vector(estimate_alpha),
                                                ppjsdm::Lightweight_square_matrix<bool>(estimate_gamma),
                                                detail::make_dispersions(model, short_range, saturation),
                                                ppjsdm::Saturated_model<double>(medium_range_model, medium_range, long_range, saturation),
                                                vector_configuration,
                                                nthreads,
                                                ppjsdm::Im_list_wrapper(covariates),
                                                initial_window_volume,
                                                restricted_window,
                                                debug);
          ++total_executions;
        }
      }
    }
    A2_plus_A3 /= total_executions;
  } else {
    detail::restrict_window(cpp_window, vector_configuration, npoints, 0.05);
    A2_plus_A3 = detail::make_A2_plus_A3(papangelou,
                                         Rcpp::as<std::vector<double>>(rho),
                                         Rcpp::as<std::vector<double>>(theta),
                                         ppjsdm::Lightweight_matrix<double>(regressors),
                                         data_list,
                                         detail::convert_list_of_boolean_matrices_to_vector(estimate_alpha),
                                         ppjsdm::Lightweight_square_matrix<bool>(estimate_gamma),
                                         detail::make_dispersions(model, short_range, saturation),
                                         ppjsdm::Saturated_model<double>(medium_range_model, medium_range, long_range, saturation),
                                         vector_configuration,
                                         nthreads,
                                         ppjsdm::Im_list_wrapper(covariates),
                                         initial_window_volume,
                                         cpp_window,
                                         debug);
  }

  return Rcpp::wrap(A2_plus_A3);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_G2_cpp(SEXP configuration,
                                   SEXP dummy,
                                   SEXP window,
                                   Rcpp::List covariates,
                                   SEXP model,
                                   Rcpp::CharacterVector medium_range_model,
                                   Rcpp::List short_range,
                                   Rcpp::NumericMatrix medium_range,
                                   Rcpp::NumericMatrix long_range,
                                   Rcpp::IntegerVector type,
                                   R_xlen_t saturation,
                                   Rcpp::NumericVector rho,
                                   Rcpp::NumericVector theta,
                                   Rcpp::NumericMatrix regressors,
                                   Rcpp::List estimate_alpha,
                                   Rcpp::LogicalMatrix estimate_gamma,
                                   int nthreads,
                                   std::string dummy_distribution,
                                   Rcpp::NumericVector mark_range) {
  // Convert the SEXP configuration to a C++ object.
  const ppjsdm::Configuration_wrapper wrapped_configuration(Rcpp::wrap(configuration));

  // Convert configuration to std::vector in order for parallelised version to work.
  const auto length_configuration(ppjsdm::size(wrapped_configuration));
  std::vector<ppjsdm::Marked_point> vector_configuration(length_configuration);
  for(decltype(ppjsdm::size(wrapped_configuration)) j(0); j < length_configuration; ++j) {
    vector_configuration[j] = wrapped_configuration[j];
  }

  // Convert the SEXP dummy to a C++ object.
  const ppjsdm::Configuration_wrapper wrapped_dummy(Rcpp::wrap(dummy));

  // Convert configuration to std::vector in order for parallelised version to work.
  const auto length_dummy(ppjsdm::size(wrapped_dummy));
  std::vector<ppjsdm::Marked_point> vector_dummy(length_dummy);
  for(decltype(ppjsdm::size(wrapped_dummy)) j(0); j < length_dummy; ++j) {
    vector_dummy[j] = wrapped_dummy[j];
  }

  // Convert the SEXP window to a C++ object.
  ppjsdm::Window cpp_window(window, mark_range);

  // Construct papangelou intensity
  const auto papangelou(detail::make_papangelou(ppjsdm::Lightweight_matrix<double>(regressors), Rcpp::as<std::vector<double>>(theta), nthreads));

  // Compute and return G2
  using g2_t = decltype(detail::make_G2_binomial(papangelou,
                                                 Rcpp::as<std::vector<double>>(rho),
                                                 ppjsdm::Lightweight_matrix<double>(regressors),
                                                 std::vector<int>{}, 1., 1));
  g2_t G2;
  if(dummy_distribution == std::string("binomial")) {
    G2 = detail::make_G2_binomial(papangelou,
                                  Rcpp::as<std::vector<double>>(rho),
                                  ppjsdm::Lightweight_matrix<double>(regressors),
                                  Rcpp::as<std::vector<int>>(type),
                                  cpp_window.volume(),
                                  nthreads);
  } else if(dummy_distribution == std::string("stratified")) {
    G2 = detail::make_G2_stratified(vector_configuration,
                                    vector_dummy,
                                    cpp_window,
                                    Rcpp::as<std::vector<double>>(theta),
                                    Rcpp::as<std::vector<double>>(rho),
                                    ppjsdm::Lightweight_matrix<double>(regressors),
                                    detail::convert_list_of_boolean_matrices_to_vector(estimate_alpha),
                                    ppjsdm::Lightweight_square_matrix<bool>(estimate_gamma),
                                    detail::make_dispersions(model, short_range, saturation),
                                    ppjsdm::Saturated_model<double>(medium_range_model, medium_range, long_range, saturation),
                                    ppjsdm::Im_list_wrapper(covariates),
                                    nthreads);
  } else {
    Rcpp::stop("Unknown dummy_distribution parameter.");
  }

  return Rcpp::wrap(G2);
}
