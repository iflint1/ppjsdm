#ifndef INCLUDE_PPJSDM_SUBSET_TO_NONNA
#define INCLUDE_PPJSDM_SUBSET_TO_NONNA

#include <Rinternals.h>

#include <algorithm> // std::remove_if
#include <type_traits> // std::remove_cv_t

namespace ppjsdm {

template<typename Configuration, typename ImList>
inline auto subset_to_nonNA(const Configuration& configuration,
                            const ImList& covariates) {
  Configuration copy(configuration);
  copy.erase(std::remove_if(copy.begin(), copy.end(),
                            [&covariates](const auto& point) {
                              for(std::remove_cv_t<decltype(covariates.size())> covariate_index(0); covariate_index < covariates.size(); ++covariate_index) {
                                const auto covariate(covariates[covariate_index](point));
                                if(R_IsNA(covariate)) {
                                  return true;
                                }
                              }
                              return false;
                            }), copy.end());

  return copy;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SUBSET_TO_NONNA
