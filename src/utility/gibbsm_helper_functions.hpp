#ifndef INCLUDE_GIBBSM_HELPERS
#define INCLUDE_GIBBSM_HELPERS

#include <Rcpp.h>

#include <string> // std::string

namespace ppjsdm {

inline auto get_number_parameters(int number_types,
                                  unsigned long long int covariates_size,
                                  Rcpp::LogicalMatrix estimate_alpha,
                                  Rcpp::LogicalMatrix estimate_gamma) {
  using size_t = unsigned long long int;
  struct Number_parameters_struct {
    size_t index_start_gamma;
    size_t index_start_covariates;
    size_t total_parameters;
  };

  size_t nalpha(0);
  size_t ngamma(0);
  for(int i(0); i < number_types; ++i) {
    for(int j(i); j < number_types; ++j) {
      if(estimate_alpha(i, j)) {
        ++nalpha;
      }
      if(estimate_gamma(i, j)) {
        ++ngamma;
      }
    }
  }

  const size_t index_start_gamma(number_types + nalpha);
  const size_t index_start_covariates(number_types + nalpha + ngamma);
  const size_t total_parameters(index_start_covariates + number_types * covariates_size);
  return Number_parameters_struct{index_start_gamma, index_start_covariates, total_parameters};
}

inline auto make_model_coloumn_names(const Im_list_wrapper& covariates,
                                     Rcpp::LogicalMatrix estimate_alpha,
                                     Rcpp::LogicalMatrix estimate_gamma) {

  const auto number_types(estimate_alpha.nrow());

  // Rcpp::CharacterVector short_range_direct_names(Rcpp::no_init(short_range_traits_input.ncol()));
  // Rcpp::CharacterVector medium_range_direct_names(Rcpp::no_init(medium_range_traits_input.ncol()));
  // Rcpp::CharacterVector short_range_joint_names(Rcpp::no_init(short_range_joint_traits_input.ncol()));
  // Rcpp::CharacterVector medium_range_joint_names(Rcpp::no_init(medium_range_joint_traits_input.ncol()));
  // for(size_t i(0); i < short_range_traits_input.ncol(); ++i) {
  //   short_range_direct_names[i] = std::string("short_range_direct_") + (i == 0 ? std::string("intercept") : std::to_string(i));
  //   medium_range_direct_names[i] = std::string("medium_range_direct_") + (i == 0 ? std::string("intercept") : std::to_string(i));
  //   short_range_joint_names[i] = std::string("short_range_joint_") + (i == 0 ? std::string("intercept") : std::to_string(i));
  //   medium_range_joint_names[i] = std::string("medium_range_joint_") + (i == 0 ? std::string("intercept") : std::to_string(i));
  // }

  //   regressors = Rcpp::no_init(log_lambda.nrow(),
  //                              log_lambda.ncol() + short_range_traits_input.ncol() + medium_range_traits_input.ncol() + short_range_joint_traits_input.ncol() + medium_range_joint_traits_input.ncol() + covariates_input.ncol());
  //   Rcpp::CharacterVector col_names(Rcpp::no_init(regressors.ncol()));
  //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(log_lambda.ncol()); ++j) {
  //     col_names[j] = log_lambda_names[j];
  //     for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
  //       regressors(i, j) = log_lambda(i, j);
  //     }
  //   }
  //   R_xlen_t index_shift(log_lambda.ncol());
  //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(short_range_traits_input.ncol()); ++j) {
  //     col_names[j + index_shift] = short_range_direct_names[j];
  //     for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
  //       regressors(i, j + index_shift) = short_range_traits_input(i, j);
  //     }
  //   }
  //   index_shift = log_lambda.ncol() + short_range_traits_input.ncol();
  //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(medium_range_traits_input.ncol()); ++j) {
  //     col_names[j + index_shift] = medium_range_direct_names[j];
  //     for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
  //       regressors(i, j + index_shift) = medium_range_traits_input(i, j);
  //     }
  //   }
  //   index_shift = log_lambda.ncol() + short_range_traits_input.ncol() + medium_range_traits_input.ncol();
  //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(short_range_joint_traits_input.ncol()); ++j) {
  //     col_names[j + index_shift] = short_range_joint_names[j];
  //     for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
  //       regressors(i, j + index_shift) = short_range_joint_traits_input(i, j);
  //     }
  //   }
  //   index_shift = log_lambda.ncol() + short_range_traits_input.ncol() + medium_range_traits_input.ncol() + short_range_joint_traits_input.ncol();
  //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(medium_range_joint_traits_input.ncol()); ++j) {
  //     col_names[j + index_shift] = medium_range_joint_names[j];
  //     for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
  //       regressors(i, j + index_shift) = medium_range_joint_traits_input(i, j);
  //     }
  //   }
  //   index_shift = log_lambda.ncol() + short_range_traits_input.ncol() + medium_range_traits_input.ncol() + short_range_joint_traits_input.ncol() + medium_range_joint_traits_input.ncol();
  //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(covariates_input.ncol()); ++j) {
  //     col_names[j + index_shift] = covariates_input_names[j];
  //     for(R_xlen_t i(0); i < regressors.nrow(); ++i) {
  //       regressors(i, j + index_shift) = covariates_input(i, j);
  //     }
  //   }
  //   Rcpp::colnames(regressors) = col_names;
  // } else {
  // Rcpp::CharacterVector col_names(Rcpp::no_init(regressors.ncol()));
  // for(R_xlen_t j(0); j < static_cast<R_xlen_t>(log_lambda.ncol()); ++j) {
  //   col_names[j] = log_lambda_names[j];
  // }
  // R_xlen_t index_shift(log_lambda.ncol());
  // if(estimate_alpha) {
  //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(alpha_input.ncol()); ++j) {
  //     col_names[j + index_shift] = alpha_names[j];
  //   }
  //   index_shift += alpha_input.ncol();
  // }
  // if(estimate_gamma) {
  //   for(R_xlen_t j(0); j < static_cast<R_xlen_t>(gamma_input.ncol()); ++j) {
  //     col_names[j + index_shift] = gamma_names[j];
  //   }
  //   index_shift += gamma_input.ncol();
  // }
  // for(R_xlen_t j(0); j < static_cast<R_xlen_t>(covariates_input.ncol()); ++j) {
  //   col_names[j + index_shift] = covariates_input_names[j];
  // }

  Rcpp::CharacterVector col_names(Rcpp::no_init(get_number_parameters(number_types,
                                                                      covariates.size(),
                                                                      estimate_alpha,
                                                                      estimate_gamma).total_parameters));
  for(R_xlen_t j(0); j < number_types; ++j) {
    col_names[j] = std::string("log_lambda") + std::to_string(j + 1);
  }
  R_xlen_t index_shift(number_types);
  for(R_xlen_t k1(0); k1 < number_types; ++k1) {
    for(R_xlen_t k2(k1); k2 < number_types; ++k2) {
      if(estimate_alpha(k1, k2)) {
        col_names[index_shift] = std::string("alpha_") + std::to_string(k1 + 1) + std::string("_") + std::to_string(k2 + 1);
        ++index_shift;
      }
    }
  }
  for(R_xlen_t k1(0); k1 < number_types; ++k1) {
    for(R_xlen_t k2(k1); k2 < number_types; ++k2) {
      if(estimate_gamma(k1, k2)) {
        col_names[index_shift] = std::string("gamma_") + std::to_string(k1 + 1) + std::string("_") + std::to_string(k2 + 1);
        ++index_shift;
      }
    }
  }
  for(R_xlen_t j(0); j < static_cast<R_xlen_t>(covariates.size()); ++j) {
    for(R_xlen_t k(0); k < number_types; ++k) {
      col_names[index_shift] = covariates.names()[j] + std::string("_") + std::to_string(k + 1);
      ++index_shift;
    }
  }

  return col_names;
}

} // namespace ppjsdm

#endif // INCLUDE_GIBBSM_HELPERS
