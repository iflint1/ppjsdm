#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration/configuration_manipulation.h"
#include "configuration/configuration_wrapper.h"
#include "configuration/get_number_points.h"
#include "phi_dispersion_model/compute_phi_dispersion.h"
#include "point/point_manipulation.h"
#include "simulation/rbinomialpp_single.h"
#include "utility/construct_if_missing.h"
#include "utility/im_wrapper.h"
#include "utility/size_t.h"
#include "utility/window_utilities.h"

#include <string> // std::string, std::to_string
#include <vector> // std::vector

template<typename Configuration>
inline Rcpp::NumericVector compute_delta_phi_dispersion(const Configuration& configuration, const ppjsdm::Marked_point& point, int number_types, Rcpp::CharacterVector model, Rcpp::NumericMatrix radius) {
  return ppjsdm::call_on_papangelou(model, radius, [&configuration, &point, number_types](const auto& papangelou) {
    return papangelou.compute(configuration, point, number_types, ppjsdm::size(configuration));
  });
}

inline void add_to_formula(std::string& formula, Rcpp::CharacterVector names) {
  for(const auto& name: names) {
    formula += std::string(" + ") + Rcpp::as<std::string>(name);
  }
}

template<typename Configuration, typename Window, typename Vector>
Rcpp::List prepare_gibbsm_data_helper(const Configuration& configuration, const Window& window, const ppjsdm::Im_list_wrapper& covariates, Rcpp::CharacterVector model, Rcpp::NumericMatrix radius, const Vector& points_by_type) {
  const auto length_configuration(ppjsdm::size(configuration));
  using size_t = ppjsdm::size_t<Configuration>;
  const size_t number_types(points_by_type.size());
  const auto covariates_length(covariates.size());
  const auto volume(window.volume());

  // Sample the dummy points D.
  // This choice of rho is the guideline from the Baddeley et al. paper, see p. 8 therein.
  // Nb: Using a lambda to inisitalise rho in order to make sure it's declared const.
  const Vector rho_times_volume([&points_by_type](){
    Vector rho_times_volume(points_by_type);
    for(auto& n: rho_times_volume) {
      n *= 4;
    }
    return rho_times_volume;
  }());
  // Set to sum(rho_times_volume), which in our case is known explicitly.
  const auto sum_rho_times_volume(4 * length_configuration);
  const auto D(ppjsdm::rbinomialpp_single<Configuration>(window, rho_times_volume, number_types, sum_rho_times_volume));
  const auto length_D(sum_rho_times_volume);

  // Default-initialise the data
  Rcpp::IntegerMatrix response(Rcpp::no_init(length_configuration + length_D, 1));
  Rcpp::NumericMatrix log_lambda(Rcpp::no_init(length_configuration + length_D, number_types));
  Rcpp::NumericMatrix rho_offset(Rcpp::no_init(length_configuration + length_D, 1));
  Rcpp::NumericMatrix alpha_input(Rcpp::no_init(length_configuration + length_D, number_types * (number_types + 1) / 2));
  Rcpp::NumericMatrix covariates_input(Rcpp::no_init(length_configuration + length_D, covariates_length * number_types));

  // Set names.
  size_t index(0);
  Rcpp::CharacterVector alpha_names(Rcpp::no_init(alpha_input.ncol()));
  for(size_t j(1); j <= number_types; ++j) {
    for(size_t k(j); k <= number_types; ++k) {
      alpha_names[index] = std::string("alpha") + std::string("_") + std::to_string(j) + std::string("_") + std::to_string(k);
      ++index;
    }
  }
  Rcpp::colnames(alpha_input) = alpha_names;
  if(covariates_length > 0) {
    Rcpp::CharacterVector covariates_input_names(Rcpp::no_init(covariates_length * number_types));
    const auto covariates_names(covariates.names());
    for(size_t i(0); i < covariates_length; ++i) {
      for(size_t j(0); j < number_types; ++j) {
        covariates_input_names[i * number_types + j] = covariates_names[i] + std::string("_") + std::to_string(j + 1);
      }
    }
    Rcpp::colnames(covariates_input) = covariates_input_names;
  }
  Rcpp::CharacterVector log_lambda_names(Rcpp::no_init(number_types));
  for(size_t i(0); i < number_types; ++i) {
    log_lambda_names[i] = std::string("log_lambda") + std::string("_") + std::to_string(i + 1);
  }
  Rcpp::colnames(log_lambda) = log_lambda_names;
  Rcpp::colnames(rho_offset) = Rcpp::wrap("rho");
  Rcpp::colnames(response) = Rcpp::wrap("response");

  // Fill the data.
  for(size_t i(0); i < length_configuration + length_D; ++i) {
    Rcpp::NumericVector dispersion;
    ppjsdm::Marked_point point;
    if(i < length_configuration) {
      response(i, 0) = 1;
      point = configuration[i];

      Configuration configuration_copy(configuration);
      ppjsdm::remove_point_by_index(configuration_copy, i);

      dispersion = compute_delta_phi_dispersion(configuration_copy, point, number_types, model, radius);
    } else {
      response(i, 0) = 0;
      point = D[i - length_configuration];
      dispersion = compute_delta_phi_dispersion(configuration, point, number_types, model, radius);
    }
    Rcpp::NumericVector location{ppjsdm::get_x(point), ppjsdm::get_y(point)};
    const size_t type_index(ppjsdm::get_type(point));

    rho_offset(i, 0) = static_cast<double>(rho_times_volume[type_index]) / volume;
    for(size_t j(0); j < covariates_length; ++j) {
      for(size_t k(0); k < number_types; ++k) {
        if(k == type_index) {
          covariates_input(i, j * number_types + k) = covariates[j](location[0], location[1]);
        } else {
          covariates_input(i, j * number_types + k) = 0;
        }
      }
    }
    R_xlen_t index(0);
    for(size_t j(0); j < number_types; ++j) {
      if(type_index == j) {
        log_lambda(i, j) = 1;
      } else {
        log_lambda(i, j) = 0;
      }
      for(size_t k(j); k < number_types; ++k) {
        if(j == type_index) {
          alpha_input(i, index) = dispersion[k];
        } else if(k == type_index) {
          alpha_input(i, index) = dispersion[j];
        } else {
          alpha_input(i, index) = 0;
        }
        ++index;
      }
    }
  }

  // Construct formula
  std::string formula("response ~ 0 + offset(-log(rho))");
  add_to_formula(formula, Rcpp::colnames(log_lambda));
  if(covariates_length > 0) {
    add_to_formula(formula, Rcpp::colnames(covariates_input));
  }
  add_to_formula(formula, Rcpp::colnames(alpha_input));

  return Rcpp::List::create(Rcpp::Named("data") = Rcpp::List::create(covariates_input,
                                        alpha_input,
                                        response,
                                        log_lambda,
                                        rho_offset),
                            Rcpp::Named("formula") = formula);
}

// [[Rcpp::export]]
Rcpp::List prepare_gibbsm_data(Rcpp::List configuration, SEXP window, Rcpp::List covariates, Rcpp::CharacterVector model, SEXP radius = R_NilValue) {
  const ppjsdm::Configuration_wrapper wrapped_configuration(configuration);
  const auto length_configuration(ppjsdm::size(wrapped_configuration));
  if(length_configuration == 0) {
    Rcpp::stop("Empty configuration.");
  }
  // The trick below allows us to find the number of different types in the configuration.
  // That number is then used to default construct `radius`.
  auto points_by_type(ppjsdm::get_number_points(wrapped_configuration));
  const auto number_types(points_by_type.size());
  radius = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(number_types, radius, 0.);
  return ppjsdm::call_on_wrapped_window(window, [&wrapped_configuration, &covariates, &model, &radius, &points_by_type](const auto& w) {
    return prepare_gibbsm_data_helper(wrapped_configuration, w, ppjsdm::Im_list_wrapper(covariates), model, radius, points_by_type);
  });
}
