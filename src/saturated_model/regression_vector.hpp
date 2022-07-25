#ifndef INCLUDE_REGRESSION_VECTOR
#define INCLUDE_REGRESSION_VECTOR

#include <vector> // std::vector

namespace ppjsdm {

template<typename ComputationType, typename Vector, typename DispersionVector, typename Dispersion, typename Matrix, typename VectorMatrices>
inline auto make_regression_vector(typename std::vector<ComputationType>::size_type type_index,
                                   typename std::vector<ComputationType>::size_type number_types,
                                   typename std::vector<ComputationType>::size_type number_parameters,
                                   typename std::vector<ComputationType>::size_type index_start_covariates,
                                   typename std::vector<ComputationType>::size_type index_start_gamma,
                                   const Vector& covariates,
                                   const DispersionVector& dispersion_short,
                                   const Dispersion& dispersion_medium,
                                   const VectorMatrices& estimate_alpha,
                                   const Matrix& estimate_gamma) {
  using size_t = typename std::vector<ComputationType>::size_type;
  std::vector<ComputationType> regression_vector(number_parameters);

  // Fill beta0
  regression_vector[type_index] = static_cast<ComputationType>(1);

  // Fill covariates
  for(decltype(covariates.size()) covariate_index(0); covariate_index < covariates.size(); ++covariate_index) {
    regression_vector[index_start_covariates + covariate_index * number_types + type_index] = static_cast<ComputationType>(covariates[covariate_index]);
  }

  // Fill alpha & gamma
  size_t index_alpha(0);
  size_t index_gamma(0);
  for(size_t j(0); j <= type_index; ++j) {
    for(size_t k(j); k < number_types; ++k) {
      for(decltype(estimate_alpha.size()) i(0); i < estimate_alpha.size(); ++i) {
        if(estimate_alpha[i](j, k)) {
          if(j == type_index) {
            regression_vector[number_types + index_alpha] = static_cast<ComputationType>(dispersion_short[i][k]);
          } else if(k == type_index) {
            regression_vector[number_types + index_alpha] = static_cast<ComputationType>(dispersion_short[i][j]);
          }
          ++index_alpha;
        }
      }
      if(estimate_gamma(j, k)) {
        if(j == type_index) {
          regression_vector[index_start_gamma + index_gamma] = static_cast<ComputationType>(dispersion_medium[k]);
        } else if(k == type_index) {
          regression_vector[index_start_gamma + index_gamma] = static_cast<ComputationType>(dispersion_medium[j]);
        }
        ++index_gamma;
      }
    }
  }

  return regression_vector;
}

} // namespace ppjsdm

#endif // INCLUDE_REGRESSION_VECTOR
