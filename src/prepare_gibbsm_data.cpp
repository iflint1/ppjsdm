#include <Rcpp.h>

#include "configuration/configuration_manipulation.h"
#include "configuration/configuration_wrapper.h"
#include "configuration/get_number_points.h"
#include "phi_dispersion_model/compute_phi_dispersion.h"
#include "point/point_manipulation.h"
#include "simulation/rbinomialpp_single.h"
#include "utility/construct_if_missing.h"
#include "utility/window_utilities.h"

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

  // This is the guideline from the Baddeley et al. paper, see p. 8 therein.
  for(auto& n: rho_times_volume) {
    n *= 4;
  }
  const auto sum_rho_times_volume(4 * length_configuration);

  const auto D(ppjsdm::rbinomialpp_single<Configuration>(window, rho_times_volume, number_types, sum_rho_times_volume));
  const auto length_D(sum_rho_times_volume);

  Rcpp::IntegerVector response(Rcpp::no_init(length_configuration + length_D));
  Rcpp::NumericMatrix log_lambda(Rcpp::no_init(length_configuration + length_D, number_types));
  Rcpp::NumericVector rho_offset(Rcpp::no_init(length_configuration + length_D));
  Rcpp::NumericMatrix alpha_input(Rcpp::no_init(length_configuration + length_D, number_types * (number_types + 1) / 2));
  Rcpp::NumericMatrix covariates_input(Rcpp::no_init(length_configuration + length_D, covariates_length));
  Rcpp::NumericMatrix traits_input(Rcpp::no_init(length_configuration + length_D, traits_length));

  for(size_t i(0); i < length_configuration + length_D; ++i) {
    // TODO: Avoid Rcpp::NumericVector
    Rcpp::NumericVector location;
    size_t type_index;
    Rcpp::NumericVector dispersion;
    if(i < length_configuration) {
      response[i] = 1;
      location = Rcpp::NumericVector{ppjsdm::get_x(configuration[i]), ppjsdm::get_y(configuration[i])};
      type_index = ppjsdm::get_type(configuration[i]);

      Configuration configuration_copy(configuration);
      ppjsdm::remove_point_by_index(configuration_copy, i);

      dispersion = compute_delta_phi_dispersion(configuration_copy, location, type_index, number_types, model, radius);
    } else {
      response[i] = 0;
      location = Rcpp::NumericVector{ppjsdm::get_x(D[i - length_configuration]), ppjsdm::get_y(D[i - length_configuration])};
      type_index = ppjsdm::get_type(D[i - length_configuration]);
      dispersion = compute_delta_phi_dispersion(configuration, location, type_index, number_types, model, radius);
    }
    rho_offset[i] = static_cast<double>(rho_times_volume[type_index]) / volume;
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
  return Rcpp::List::create(Rcpp::Named("response") = response,
                            Rcpp::Named("log_lambda") = log_lambda,
                            Rcpp::Named("rho_offset") = rho_offset,
                            Rcpp::Named("alpha_input") = alpha_input,
                            Rcpp::Named("covariates_input") = covariates_input,
                            Rcpp::Named("traits_input") = traits_input);
  // Rcpp::List additional_traits;
  //
  //
  // if(ncovariates > 0) {
  //   if(ntraits > 0) {
  //     additional_traits = Rcpp::append(covariate_list, trait_list);
  //   } else {
  //     additional_traits = Rcpp::append(covariate_list, alpha_list);
  //   }
  //
  // } else {
  //   if(ntraits > 0) {
  //     additional_traits = trait_list;
  //   } else {
  //     additional_traits = alpha_list;
  //   }
  // }
  // Rcpp::DataFrame data(Rcpp::append(additional_traits, Rcpp::List(Rcpp::Named("response") = response, Rcpp::Named("log_lambda") = log_lambda, Rcpp::Named("rho") = rho_offset)));
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
