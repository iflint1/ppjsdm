#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_manipulation.hpp"
#include "configuration/configuration_wrapper.hpp"
#include "configuration/get_number_points.hpp"

#include "point/point_manipulation.hpp"

#include "saturated_varphi_model/saturated_varphi_model.hpp"

#include "simulation/rbinomialpp_single.hpp"

#include "utility/construct_if_missing.hpp"
#include "utility/im_wrapper.hpp"
#include "utility/lightweight_matrix.hpp"
#include "utility/window_utilities.hpp"

#include <algorithm> // std::remove_if
#include <cmath> // std::log
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

template<typename Configuration, typename Window, typename DispersionModel, typename Vector>
Rcpp::List prepare_gibbsm_data_helper(const Configuration& configuration, const Window& window, const ppjsdm::Im_list_wrapper& covariates, const DispersionModel& dispersion_model, const Vector& points_by_type) {
  const auto length_configuration(ppjsdm::size(configuration));
  using size_t = ppjsdm::size_t<Configuration>;
  const size_t number_types(points_by_type.size());
  const size_t covariates_length(covariates.size());
  const auto volume(window.volume());

  // Sample the dummy points D.
  // This choice of rho is the guideline from the Baddeley et al. paper, see p. 8 therein.
  Vector rho_times_volume(points_by_type);
  size_t length_rho_times_volume(0);
  for(auto& n: rho_times_volume) {
    const auto mult_by_four(n * 4);
    n = mult_by_four < 500 ? 500 : mult_by_four;
    length_rho_times_volume += n;
  }
  auto D(ppjsdm::rbinomialpp_single<std::vector<ppjsdm::Marked_point>>(window, rho_times_volume, number_types, length_rho_times_volume));


  // Get rid of locations with an NA value for one of the covariates
  // TODO: I'm recomputing all the covariates below; might be able to remove the points at that point!
  D.erase(std::remove_if(D.begin(), D.end(), [&covariates, covariates_length](const auto& point) {
    for(size_t j(0); j < covariates_length; ++j) {
      const auto covariate(covariates[j](ppjsdm::get_x(point),  ppjsdm::get_y(point)));
      if(R_IsNA(covariate)) {
        return true;
      }
    }
    return false;
  }), D.end());
  const size_t length_D(D.size());

  // Default-initialise the data
  std::vector<int> response(length_configuration + length_D);
  std::vector<double> rho_offset(length_configuration + length_D);
  ppjsdm::Lightweight_matrix<double> log_lambda(length_configuration + length_D, number_types);
  ppjsdm::Lightweight_matrix<double> alpha_input(length_configuration + length_D, number_types * (number_types + 1) / 2);
  ppjsdm::Lightweight_matrix<double> covariates_input(length_configuration + length_D, covariates_length * number_types);

  // Make shift vector
  Rcpp::NumericVector shift(Rcpp::no_init(number_types));
  for(size_t j(0); j < number_types; ++j) {
    shift[j] = -std::log(static_cast<double>(rho_times_volume[j]) / volume);
  }

  // Precompute dispersion and other computation-intensive things.
  struct Result_type {
    Result_type(): is_in_configuration(), type(), dispersion(), covariates() {}

    Result_type(bool i, int t, std::vector<double> v, std::vector<double> w):
      is_in_configuration(i), type(t), dispersion(std::move(v)), covariates(std::move(w)) {}

    bool is_in_configuration;
    int type;
    std::vector<double> dispersion;
    std::vector<double> covariates;
  };
  std::vector<Result_type> precomputed_results;
  precomputed_results.reserve(length_configuration + length_D);

#pragma omp parallel
{
  std::vector<Result_type> results_private;
  results_private.reserve(length_configuration + length_D);
#pragma omp for nowait
  for(size_t i = 0; i < length_configuration; ++i) {
    // TODO: Avoidable?
    Configuration configuration_copy(configuration);
    ppjsdm::remove_point_by_index(configuration_copy, i);
    const auto d(dispersion_model.compute(configuration[i], number_types, configuration_copy));
    std::vector<double> cov(covariates_length);
    for(size_t k(0); k < covariates_length; ++k) {
      const auto covariate(covariates[k](configuration[i]));
      if(R_IsNA(covariate)) {
        Rcpp::stop("One of the covariates' value is NA on one of the locations in the dataset.");
      }
      cov[k] = covariate;
    }
    results_private.emplace_back(true, ppjsdm::get_type(configuration[i]), std::move(d), std::move(cov));
  }
#pragma omp for nowait
  for(size_t i = 0; i < length_D; ++i) {
    const auto d(dispersion_model.compute(D[i], number_types, configuration));
    std::vector<double> cov(covariates_length);
    for(size_t k(0); k < covariates_length; ++k) {
      const auto covariate(covariates[k](D[i]));
      if(R_IsNA(covariate)) {
        Rcpp::stop("One of the covariates' value is NA on one of the locations in the dataset.");
      }
      cov[k] = covariate;
    }
    results_private.emplace_back(false, ppjsdm::get_type(D[i]), std::move(d), std::move(cov));
  }
#pragma omp critical
  precomputed_results.insert(precomputed_results.end(), results_private.begin(), results_private.end());
}

  // Fill the regressors, response, offset and shift with what we precomputed.
  for(size_t i(0); i < precomputed_results.size(); ++i) {
    response[i] = precomputed_results[i].is_in_configuration ? 1 : 0;
    const size_t type_index(precomputed_results[i].type);

    rho_offset[i] = -std::log(static_cast<double>(rho_times_volume[type_index]) / volume);
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
          alpha_input(i, index++) = precomputed_results[i].dispersion[k];
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
            alpha_input(i, index++) = precomputed_results[i].dispersion[j];
          } else {
            alpha_input(i, index++) = 0;
          }
        }
      }
    }
  }

  // Convert response and rho_offset to Rcpp objects.
  Rcpp::IntegerMatrix response_rcpp(Rcpp::no_init(response.size(), 1));
  for(R_xlen_t i(0); i < static_cast<R_xlen_t>(response.size()); ++i) {
    response_rcpp(i, 0) = response[i];
  }
  Rcpp::NumericMatrix rho_offset_rcpp(Rcpp::no_init(rho_offset.size(), 1));
  for(R_xlen_t i(0); i < static_cast<R_xlen_t>(rho_offset.size()); ++i) {
    rho_offset_rcpp(i, 0) = rho_offset[i];
  }

  // Set names.
  Rcpp::CharacterVector alpha_names(Rcpp::no_init(alpha_input.ncol()));
  size_t index(0);
  for(size_t j(0); j < number_types; ++j) {
    for(size_t k(j); k < number_types; ++k) {
      alpha_names[index++] = std::string("alpha_") + std::to_string(j + 1) + std::string("_") + std::to_string(k + 1);
    }
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

  // TODO: Write rbind that also works with names?
  // TODO: Might also be simpler to just construct regressors straight away
  // Put all the regressors into a unique matrix that'll be sent to glm / glmnet.
  Rcpp::NumericMatrix regressors(Rcpp::no_init(log_lambda.nrow(),
                                               log_lambda.ncol() + alpha_input.ncol() + covariates_input.ncol()));
  Rcpp::CharacterVector col_names(Rcpp::no_init(regressors.ncol()));
  for(R_xlen_t j(0); j < static_cast<R_xlen_t>(log_lambda.ncol()); ++j) {
    col_names[j] = log_lambda_names[j];
    for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
      regressors(i, j) = log_lambda(i, j);
    }
  }
  for(R_xlen_t j(0); j < static_cast<R_xlen_t>(alpha_input.ncol()); ++j) {
    col_names[j + log_lambda.ncol()] = alpha_names[j];
    for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
      regressors(i, j + log_lambda.ncol()) = alpha_input(i, j);
    }
  }
  for(R_xlen_t j(0); j < static_cast<R_xlen_t>(covariates_input.ncol()); ++j) {
    col_names[j + log_lambda.ncol() + alpha_input.ncol()] = covariates_input_names[j];
    for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
      regressors(i, j + log_lambda.ncol() + alpha_input.ncol()) = covariates_input(i, j);
    }
  }
  Rcpp::colnames(regressors) = col_names;


  return Rcpp::List::create(Rcpp::Named("response") = response_rcpp,
                            Rcpp::Named("offset") = rho_offset_rcpp,
                            Rcpp::Named("regressors") = regressors,
                            Rcpp::Named("shift") = shift
                            );
}

// [[Rcpp::export]]
Rcpp::List prepare_gibbsm_data(SEXP configuration, SEXP window, Rcpp::List covariates, Rcpp::CharacterVector model, Rcpp::CharacterVector medium_range_model, SEXP short_range, SEXP medium_range, SEXP long_range, R_xlen_t saturation) {
  const ppjsdm::Configuration_wrapper wrapped_configuration(configuration);
  const auto length_configuration(ppjsdm::size(wrapped_configuration));
  if(length_configuration == 0) {
    Rcpp::stop("Empty configuration.");
  }
  // Convert to std::vector in order for parallelised version to work.
  std::vector<ppjsdm::Marked_point> vector_configuration(length_configuration);
  for(decltype(ppjsdm::size(wrapped_configuration)) i(0); i < length_configuration; ++i) {
    vector_configuration[i] = wrapped_configuration[i];
  }
  return ppjsdm::call_on_wrapped_window(window, [&vector_configuration, &covariates, model, medium_range_model, short_range, medium_range, long_range, saturation](const auto& w) {
    // The trick below allows us to find the number of different types in the configuration.
    // That number is then used to default construct `short_range`.
    const auto points_by_type(ppjsdm::get_number_points(vector_configuration));
    const auto number_types(points_by_type.size());
    const auto sh(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(short_range, 0.1 * w.diameter(), number_types));
    return ppjsdm::call_on_dispersion_model(model, sh, saturation, [number_types, medium_range_model, medium_range, long_range, saturation, &vector_configuration, &w, &covariates, &points_by_type](const auto& short_papangelou) {
      // TODO: Use medium range Papangelou
      const auto me(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(medium_range, 0.1 * w.diameter(), number_types));
      const auto lo(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(long_range, 0.2 * w.diameter(), number_types));
      return ppjsdm::call_on_medium_range_dispersion_model(medium_range_model, me, lo, saturation, [&short_papangelou, &vector_configuration, &w, &covariates, &points_by_type](const auto&) {
        return prepare_gibbsm_data_helper(vector_configuration, w, ppjsdm::Im_list_wrapper(covariates), short_papangelou, points_by_type);
      });
    });
  });
}
