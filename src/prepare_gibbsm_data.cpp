#include <Rcpp.h>

#include "configuration/configuration_manipulation.h"
#include "configuration/configuration_wrapper.h"
#include "configuration/get_number_points.h"
#include "phi_dispersion_model/compute_phi_dispersion.h"
#include "point/point_manipulation.h"
#include "simulation/rbinomialpp_single.h"
#include "utility/construct_if_missing.h"
#include "utility/window_utilities.h"

#include <string> // std::string, std::to_string
#include <type_traits> // std::remove_const
#include <vector> // std::vector

// TODO: Use marked point
template<typename Configuration>
Rcpp::NumericVector compute_delta_phi_dispersion(const Configuration& configuration, Rcpp::NumericVector location, R_xlen_t type, int number_types, Rcpp::CharacterVector model = "identity", Rcpp::NumericMatrix radius = Rcpp::NumericMatrix(1, 1)) {
  return ppjsdm::call_on_papangelou(model, radius, [&configuration, &location, type, number_types](const auto& papangelou) {
    const auto number_points(ppjsdm::size(configuration));
    return papangelou.compute(configuration, ppjsdm::Marked_point(location[0], location[1], type), number_types, number_points);
  });
}

template<typename Configuration, typename Window>
Rcpp::List prepare_gibbsm_data_helper(const Configuration& configuration, const Window& window, Rcpp::List covariates, Rcpp::List traits, Rcpp::CharacterVector model, Rcpp::NumericMatrix radius) {
  const auto length_configuration(ppjsdm::size(configuration));
  if(length_configuration == 0) {
    Rcpp::stop("Empty configuration.");
  }
  auto rho_times_volume(ppjsdm::get_number_points(configuration));
  const auto number_types(rho_times_volume.size());
  using size_t = std::remove_const_t<decltype(number_types)>;
  const auto covariates_length(covariates.size());
  const auto traits_length(traits.size());
  const auto volume(window.volume());

  // Sample the dummy points D.
  // This choice of rho is the guideline from the Baddeley et al. paper, see p. 8 therein.
  for(auto& n: rho_times_volume) {
    n *= 4;
  }
  const auto sum_rho_times_volume(4 * length_configuration);
  const auto D(ppjsdm::rbinomialpp_single<Configuration>(window, rho_times_volume, number_types, sum_rho_times_volume));
  const auto length_D(sum_rho_times_volume);

  // Default-initialise the data
  Rcpp::IntegerMatrix response(Rcpp::no_init(length_configuration + length_D, 1));
  Rcpp::NumericMatrix log_lambda(Rcpp::no_init(length_configuration + length_D, number_types));
  Rcpp::NumericMatrix rho_offset(Rcpp::no_init(length_configuration + length_D, 1));
  Rcpp::NumericMatrix alpha_input(Rcpp::no_init(length_configuration + length_D, number_types * (number_types + 1) / 2));
  Rcpp::NumericMatrix covariates_input(Rcpp::no_init(length_configuration + length_D, covariates_length));
  Rcpp::NumericMatrix traits_input(Rcpp::no_init(length_configuration + length_D, traits_length));

  // Set names.
  size_t index(0);
  Rcpp::CharacterVector alpha_names(Rcpp::no_init(alpha_input.ncol()));
  for(size_t j(1); j <= number_types; ++j) {
    for(size_t k(j); k <= number_types; ++k) {
      alpha_names[index] = std::string("alpha") + std::to_string(j) + std::string("_") + std::to_string(k);
      ++index;
    }
  }
  Rcpp::colnames(alpha_input) = alpha_names;
  if(covariates_length > 0) {
    Rcpp::colnames(covariates_input) = Rcpp::as<Rcpp::CharacterVector>(covariates.names());
  }
  if(traits_length > 0) {
    Rcpp::colnames(traits_input) = Rcpp::as<Rcpp::CharacterVector>(traits.names());
  }
  Rcpp::CharacterVector log_lambda_names(Rcpp::no_init(number_types));
  for(size_t i(0); i < number_types; ++i) {
    log_lambda_names[i] = std::string("log_lambda") + std::to_string(i + 1);
  }
  Rcpp::colnames(log_lambda) = log_lambda_names;
  Rcpp::colnames(rho_offset) = Rcpp::wrap("rho");
  Rcpp::colnames(response) = Rcpp::wrap("response");

  // Fill the data.
  for(size_t i(0); i < length_configuration + length_D; ++i) {
    // TODO: Avoid Rcpp::NumericVector
    Rcpp::NumericVector location;
    size_t type_index;
    Rcpp::NumericVector dispersion;
    if(i < length_configuration) {
      response(i, 0) = 1;
      location = Rcpp::NumericVector{ppjsdm::get_x(configuration[i]), ppjsdm::get_y(configuration[i])};
      type_index = ppjsdm::get_type(configuration[i]);

      Configuration configuration_copy(configuration);
      ppjsdm::remove_point_by_index(configuration_copy, i);

      dispersion = compute_delta_phi_dispersion(configuration_copy, location, type_index, number_types, model, radius);
    } else {
      response(i, 0) = 0;
      location = Rcpp::NumericVector{ppjsdm::get_x(D[i - length_configuration]), ppjsdm::get_y(D[i - length_configuration])};
      type_index = ppjsdm::get_type(D[i - length_configuration]);
      dispersion = compute_delta_phi_dispersion(configuration, location, type_index, number_types, model, radius);
    }
    rho_offset(i, 0) = static_cast<double>(rho_times_volume[type_index]) / volume;
    for(R_xlen_t j(0); j < covariates_length; ++j) {
      // TODO: Use im access function
      covariates_input(i, j) = Rf_asReal(Rcpp::as<Rcpp::Function>(covariates[j])(location[0], location[1]));
    }
    for(R_xlen_t j(0); j < traits_length; ++j) {
      double inner_product(0);
      for(size_t k(0); k < number_types; ++k) {
        inner_product += dispersion[k] * Rcpp::as<Rcpp::NumericMatrix>(traits[j])(type_index, k);
      }
      traits_input(i, j) = inner_product;
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

  std::string formula("response ~ 0 + offset(-log(rho))");

  const auto log_lambda_col_names = Rcpp::as<Rcpp::CharacterVector>(Rcpp::colnames(log_lambda));
  for(R_xlen_t i(0); i < log_lambda_col_names.size(); ++i) {
    formula += std::string(" + ") + Rcpp::as<std::string>(log_lambda_col_names[i]);
  }

  if(covariates_length > 0) {
    const auto covariates_input_col_names = Rcpp::as<Rcpp::CharacterVector>(Rcpp::colnames(covariates_input));
    for(R_xlen_t i(0); i < covariates_input_col_names.size(); ++i) {
      formula += std::string(" + ") + Rcpp::as<std::string>(covariates_input_col_names[i]);
    }
  }

  if(traits_length > 0) {
    const auto traits_input_col_names = Rcpp::as<Rcpp::CharacterVector>(Rcpp::colnames(traits_input));
    for(R_xlen_t i(0); i < traits_input_col_names.size(); ++i) {
      formula += std::string(" + ") + Rcpp::as<std::string>(traits_input_col_names[i]);
    }
  } else {
    const auto alpha_col_names = Rcpp::as<Rcpp::CharacterVector>(Rcpp::colnames(alpha_input));
    for(R_xlen_t i(0); i < alpha_col_names.size(); ++i) {
      formula += std::string(" + ") + Rcpp::as<std::string>(alpha_col_names[i]);
    }
  }

  return Rcpp::List::create(Rcpp::Named("response") = response,
                            Rcpp::Named("log_lambda") = log_lambda,
                            Rcpp::Named("rho") = rho_offset,
                            Rcpp::Named("alpha") = alpha_input,
                            Rcpp::Named("covariates") = covariates_input,
                            Rcpp::Named("traits") = traits_input,
                            Rcpp::Named("formula") = formula);
}

// [[Rcpp::export]]
Rcpp::List prepare_gibbsm_data(Rcpp::List configuration, SEXP window = R_NilValue, Rcpp::List covariates = Rcpp::List(), Rcpp::List traits = Rcpp::List(), Rcpp::CharacterVector model = "identity", SEXP radius = R_NilValue) {
  const ppjsdm::Configuration_wrapper wrapped_configuration(configuration);
  const auto length_configuration(ppjsdm::size(wrapped_configuration));
  if(length_configuration == 0) {
    Rcpp::stop("Empty configuration.");
  }
  std::vector<ppjsdm::Marked_point> vector_configuration(length_configuration);

  for(R_xlen_t i(0); i < length_configuration; ++i) {
    vector_configuration[i] = wrapped_configuration[i];
  }
  // TODO: Remove this part
  auto rho_times_volume(ppjsdm::get_number_points(vector_configuration));
  const auto number_types(rho_times_volume.size());
  radius = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(number_types, radius, 0.);
  return ppjsdm::call_on_wrapped_window(window, [&vector_configuration, &covariates, &traits, &model, &radius](const auto& w) {
    return prepare_gibbsm_data_helper(vector_configuration, w, covariates, traits, model, radius);
  });
}
