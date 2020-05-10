#include <Rcpp.h>
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

// TODO: I suspect max_points_by_type isn't needed, only need number_types
template<typename Configuration, typename Vector, typename Model>
Rcpp::List compute_A2_A3_helper(const Configuration& configuration, const ppjsdm::Im_list_wrapper& covariates, const ppjsdm::Saturated_model& dispersion_model, const ppjsdm::Saturated_model& medium_dispersion_model, const Model& model, const Vector& max_points_by_type, double rho, double area, Rcpp::List t_over_papangelou) {
  using size_t = ppjsdm::size_t<Configuration>;
  const int number_types(max_points_by_type.size());

  const size_t total_points(ppjsdm::size(configuration));
  const size_t total_parameters(number_types * (2 + covariates.size() + number_types));

  Rcpp::NumericMatrix A2(total_parameters, total_parameters);
  Rcpp::NumericMatrix A3(total_parameters, total_parameters);

  for(size_t i(0); i < total_points; ++i) {
    const int type_i(ppjsdm::get_type(configuration[i]));

    std::vector<double> cov_i(covariates.size());
    for(R_xlen_t k(0); k < covariates.size(); ++k) {
      cov_i[k] = covariates[k](configuration[i]);
    }

    Rcpp::NumericVector t_over_papangelou_i(0);
    for(R_xlen_t l(0); l < t_over_papangelou.size(); ++l) {
      const auto current(Rcpp::as<Rcpp::List>(t_over_papangelou[l]));
      if(ppjsdm::is_equal(configuration[i], ppjsdm::Marked_point(current["x"], current["y"], Rcpp::as<int>(current["type"]) - 1, current["mark"]))) {
        t_over_papangelou_i = Rcpp::as<Rcpp::NumericVector>(current["value"]);
        break;
      }
    }
    if(t_over_papangelou_i.size() == 0) {
      Rcpp::stop("Did not find the current point in t_over_papangelou");
    }

    for(size_t j(0); j < total_points; ++j) {
      if(j != i) {
        const int type_j(ppjsdm::get_type(configuration[j]));

        Rcpp::NumericVector t_over_papangelou_j(0);
        for(R_xlen_t l(0); l < t_over_papangelou.size(); ++l) {
          const auto current(Rcpp::as<Rcpp::List>(t_over_papangelou[l]));
          if(ppjsdm::is_equal(configuration[j], ppjsdm::Marked_point(current["x"], current["y"], Rcpp::as<int>(current["type"]) - 1, current["mark"]))) {
            t_over_papangelou_j = Rcpp::as<Rcpp::NumericVector>(current["value"]);
            break;
          }
        }
        if(t_over_papangelou_j.size() == 0) {
          Rcpp::stop("Did not find the current point in t_over_papangelou");
        }

        Configuration configuration_copy(configuration);
        ppjsdm::remove_point(configuration_copy, configuration[j]);
        const auto papangelou_j_minus_j(model.compute_papangelou(configuration[j], configuration_copy));
        ppjsdm::remove_point(configuration_copy, configuration[i]);
        const auto papangelou_i_minus_two(model.compute_papangelou(configuration[i], configuration_copy));
        const auto papangelou_j_minus_two(model.compute_papangelou(configuration[j], configuration_copy));

        const auto short_i(ppjsdm::compute_dispersion(dispersion_model, configuration[i], number_types, configuration_copy));
        const auto medium_i(ppjsdm::compute_dispersion(medium_dispersion_model, configuration[i], number_types, configuration_copy));

        const auto short_j(ppjsdm::compute_dispersion(dispersion_model, configuration[j], number_types, configuration_copy));
        const auto medium_j(ppjsdm::compute_dispersion(medium_dispersion_model, configuration[j], number_types, configuration_copy));

        std::vector<double> cov_j(covariates.size());
        for(R_xlen_t k(0); k < covariates.size(); ++k) {
          cov_j[k] = covariates[k](configuration[j]);
        }

        std::vector<double> t_i(total_parameters);
        std::vector<double> t_j(total_parameters);

        size_t index(number_types);
        for(int k1(0); k1 < number_types; ++k1) {
          if(k1 == type_i) {
            t_i[k1] = 1.;

            for(int k2(0); k2 < covariates.size(); ++k2) {
              t_i[number_types * (2 + number_types) + k2 * number_types + k1] = cov_i[k2];
            }

            for(int k2(k1); k2 < number_types; ++k2) {
              t_i[index] = short_i[k2];
              t_i[index + number_types * (number_types + 1) / 2] = medium_i[k2];
            }
          } else {
            t_i[k1] = 0.;

            for(int k2(0); k2 < covariates.size(); ++k2) {
              t_i[number_types * (2 + number_types) + k2 * number_types + k1] = 0.;
            }

            for(int k2(k1); k2 < number_types; ++k2) {
              if(k2 == type_i) {
                t_i[index] = short_i[k1];
                t_i[index + number_types * (number_types + 1) / 2] = medium_i[k1];
              } else {
                t_i[index] = 0.;
                t_i[index + number_types * (number_types + 1) / 2] = 0.;
              }
            }
          }

          if(k1 == type_j) {
            t_j[k1] = 1.;

            for(int k2(0); k2 < covariates.size(); ++k2) {
              t_j[number_types * (2 + number_types) + k2 * number_types + k1] = cov_j[k2];
            }

            for(int k2(k1); k2 < number_types; ++k2) {
              t_j[index] = short_j[k2];
              t_j[index + number_types * (number_types + 1) / 2] = medium_j[k2];
              ++index;
            }
          } else {
            t_j[k1] = 0.;

            for(int k2(0); k2 < covariates.size(); ++k2) {
              t_j[number_types * (2 + number_types) + k2 * number_types + k1] = 0.;
            }

            for(int k2(k1); k2 < number_types; ++k2) {
              if(k2 == type_j) {
                t_j[index] = short_j[k1];
                t_j[index + number_types * (number_types + 1) / 2] = medium_j[k1];
              } else {
                t_j[index] = 0.;
                t_j[index + number_types * (number_types + 1) / 2] = 0.;
              }
              ++index;
            }
          }
        }

        for(size_t k1(0); k1 < total_parameters; ++k1) {
          for(size_t k2(0); k2 < total_parameters; ++k2) {
            const auto result(t_i[k1] * t_j[k2] / ((papangelou_i_minus_two + rho) * (papangelou_j_minus_two + rho)) * (papangelou_j_minus_two / papangelou_j_minus_j - 1.));
            A2(k1, k2) += result;
            A3(k1, k2) += (t_over_papangelou_i[k1] - t_i[k1] / (papangelou_i_minus_two + rho)) * (t_over_papangelou_j[k2] - t_j[k2] / (papangelou_j_minus_two + rho));
          }
        }
      }

      // TODO: index or formula here and in other vectors?
      /*size_t index(0);
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
            alpha_input(i, index) = precomputed_results[i].dispersion[k];
            gamma_input(i, index++) = precomputed_results[i].medium_dispersion[k];
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
              alpha_input(i, index) = precomputed_results[i].dispersion[j];
              gamma_input(i, index++) = precomputed_results[i].medium_dispersion[j];
            } else {
              // TODO: Not necessary since zero-initialized.
              alpha_input(i, index) = 0;
              gamma_input(i, index++) = 0;
            }
          }
        }
      }*/
    }
  }

  for(size_t k1(0); k1 < total_parameters; ++k1) {
    for(size_t k2(0); k2 < total_parameters; ++k2) {
      A2(k1, k2) *= rho * rho / area;
      A3(k1, k2) *= rho * rho / area;
    }
  }

  // Set names
  Rcpp::CharacterVector col_names(Rcpp::no_init(A2.ncol()));
  for(R_xlen_t j(0); j < number_types; ++j) {
    col_names[j] = std::string("log_lambda") + std::to_string(j + 1);
  }
  R_xlen_t index_shift(number_types);
  for(R_xlen_t k1(0); k1 < number_types; ++k1) {
    for(R_xlen_t k2(k1); k2 < number_types; ++k2) {
      col_names[index_shift] = std::string("alpha_") + std::to_string(k1 + 1) + std::string("_") + std::to_string(k2 + 1);
      ++index_shift;
    }
  }
  for(R_xlen_t k1(0); k1 < number_types; ++k1) {
    for(R_xlen_t k2(k1); k2 < number_types; ++k2) {
      col_names[index_shift] = std::string("gamma_") + std::to_string(k1 + 1) + std::string("_") + std::to_string(k2 + 1);
      ++index_shift;
    }
  }
  for(R_xlen_t j(0); j < static_cast<R_xlen_t>(covariates.size()); ++j) {
    col_names[index_shift] = covariates.names()[j];
    ++index_shift;
  }
  Rcpp::colnames(A2) = col_names;
  Rcpp::rownames(A2) = col_names;

  Rcpp::colnames(A3) = col_names;
  Rcpp::rownames(A3) = col_names;

  return Rcpp::List::create(Rcpp::Named("A2") = A2,
                            Rcpp::Named("A3") = A3);
}

// [[Rcpp::export]]
Rcpp::List compute_A2_A3(SEXP configuration, Rcpp::List covariates, Rcpp::CharacterVector model, Rcpp::CharacterVector medium_range_model, SEXP short_range, SEXP medium_range, SEXP long_range, R_xlen_t saturation, Rcpp::NumericMatrix alpha, Rcpp::NumericVector lambda, Rcpp::NumericMatrix beta, Rcpp::NumericMatrix gamma, double rho, double area, Rcpp::List t_over_papangelou) {
  // Construct std::vector of configurations.
  const ppjsdm::Configuration_wrapper wrapped_configuration(Rcpp::wrap(configuration));
  const auto length_configuration(ppjsdm::size(wrapped_configuration));

  // Convert configurations to std::vector in order for parallelised version to work.
  std::vector<ppjsdm::Marked_point> vector_configuration(length_configuration);
  for(decltype(ppjsdm::size(wrapped_configuration)) j(0); j < length_configuration; ++j) {
    vector_configuration[j] = wrapped_configuration[j];
  }

  // The trick below allows us to find the number of different types in the configuration.
  // That number is then used to default construct `short_range`.
  const auto points_by_type(ppjsdm::get_number_points(vector_configuration));

  const auto dispersion(ppjsdm::Saturated_model(model, short_range, saturation));
  const auto medium_range_dispersion(ppjsdm::Saturated_model(medium_range_model, medium_range, long_range, saturation));

  const ppjsdm::Truncated_exponential_family_model<Rcpp::NumericVector> exponential_model(lambda,
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

  return compute_A2_A3_helper(vector_configuration, ppjsdm::Im_list_wrapper(covariates), dispersion, medium_range_dispersion, exponential_model, points_by_type, rho, area, t_over_papangelou);
}
