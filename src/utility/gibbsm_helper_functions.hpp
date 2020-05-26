#ifndef INCLUDE_GIBBSM_HELPERS
#define INCLUDE_GIBBSM_HELPERS

#include <Rcpp.h>

#include <string> // std::string

namespace ppjsdm {

inline auto get_number_parameters(int number_types,
                                  unsigned long long int covariates_size,
                                  bool estimate_alpha,
                                  bool estimate_gamma) {
  using size_t = unsigned long long int;
  struct Number_parameters_struct {
    size_t index_start_gamma;
    size_t index_start_covariates;
    size_t total_parameters;
  };

  size_t index_start_gamma(0);
  size_t index_start_covariates(0);
  if(estimate_alpha) {
    if(estimate_gamma) {
      index_start_gamma = number_types + number_types * (number_types + 1) / 2;
      index_start_covariates = number_types * (2 + number_types);
    } else {
      index_start_covariates = number_types + number_types * (number_types + 1) / 2;
    }
  } else {
    if(estimate_gamma) {
      index_start_gamma = number_types;
      index_start_covariates = number_types + number_types * (number_types + 1) / 2;
    } else {
      index_start_covariates = number_types;
    }
  }
  const size_t total_parameters(index_start_covariates + number_types * covariates_size);
  return Number_parameters_struct{index_start_gamma, index_start_covariates, total_parameters};
}

inline auto make_model_coloumn_names(const Im_list_wrapper& covariates,
                                     int number_types,
                                     bool estimate_alpha,
                                     bool estimate_gamma) {
  Rcpp::CharacterVector col_names(Rcpp::no_init(get_number_parameters(number_types, covariates.size(), estimate_alpha, estimate_gamma).total_parameters));
  for(R_xlen_t j(0); j < number_types; ++j) {
    col_names[j] = std::string("log_lambda") + std::to_string(j + 1);
  }
  R_xlen_t index_shift(number_types);
  if(estimate_alpha) {
    for(R_xlen_t k1(0); k1 < number_types; ++k1) {
      for(R_xlen_t k2(k1); k2 < number_types; ++k2) {
        col_names[index_shift] = std::string("alpha_") + std::to_string(k1 + 1) + std::string("_") + std::to_string(k2 + 1);
        ++index_shift;
      }
    }
  }
  if(estimate_gamma) {
    for(R_xlen_t k1(0); k1 < number_types; ++k1) {
      for(R_xlen_t k2(k1); k2 < number_types; ++k2) {
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
