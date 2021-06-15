#ifndef INCLUDE_COMPUTE_DISPERSION_FITTING
#define INCLUDE_COMPUTE_DISPERSION_FITTING

#include <Rcpp.h>

#include "compute_dispersion_implementation.hpp"
#include "saturated_model.hpp"

#include <type_traits> // std::remove_cv_t
#include <utility> // std::forward
#include <vector> // std::vector

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ppjsdm {
namespace detail {

template<typename AbstractDispersion>
struct dispersion_computation_fitting {
  template<typename Configuration, typename FloatType, typename OtherConfiguration>
  auto operator()(const Saturated_model<FloatType>& varphi,
                R_xlen_t number_types,
                int nthreads,
                bool compute_on_configuration,
                const Configuration& configuration,
                const OtherConfiguration& other_configuration) const {
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif

    using ValueType = typename AbstractDispersion::ValueType;
    using CountType = std::vector<ValueType>;
    using DispersionType = std::vector<FloatType>;

    const auto configuration_size(size(configuration));
    using size_t = std::remove_cv_t<decltype(size(configuration))>;
    std::vector<CountType> count_vector(configuration_size);
    const auto other_configuration_size(size(other_configuration));
    const auto index_configuration(compute_on_configuration ? configuration_size : 0);
    std::vector<DispersionType> dispersion(index_configuration + other_configuration_size);

    // TODO: The code below seems to work, but it can perhaps be optimized or grouped together with the non-saturated
    // case below.
    if(static_cast<decltype(size(configuration))>(varphi.get_saturation()) >= size(configuration)) {
      for(size_t i(0); i < configuration_size; ++i) {
        count_vector[i] = CountType(number_types);
        for(size_t j(0); j < i; ++j) {
          // TODO: varphi(configuration[i], configuration[j]) only needs to be computed once
          AbstractDispersion::template update_count<std::numeric_limits<int>::max()>(varphi, count_vector[i][get_type(configuration[j])],
                                                                                  configuration[i], configuration[j]);
          AbstractDispersion::template update_count<std::numeric_limits<int>::max()>(varphi, count_vector[j][get_type(configuration[i])],
                                                                                  configuration[i], configuration[j]);
        }
      }

      if(compute_on_configuration) {
        for(size_t i(0); i < configuration_size; ++i) {
          dispersion[i] = DispersionType(number_types);
          add_count_to_dispersion<0, AbstractDispersion, 2>(varphi, dispersion[i], count_vector[i], configuration[i]);
        }
      }

#pragma omp parallel default(none) shared(other_configuration) \
      shared(dispersion, configuration, number_types, count_vector, varphi)
{
      decltype(dispersion) dispersion_private(size(other_configuration));
#pragma omp for nowait
      for(decltype(size(other_configuration)) i = 0; i < size(other_configuration); ++i) {
        dispersion_private[i] = DispersionType(number_types);
        CountType count_point(number_types);
        for(size_t j(0); j < configuration_size; ++j) {
          AbstractDispersion::template update_count<std::numeric_limits<int>::max()>(varphi, count_point[get_type(configuration[j])],
                                                                                  other_configuration[i], configuration[j]);
        }
        add_count_to_dispersion<0, AbstractDispersion, 2>(varphi, dispersion_private[i],
                                                          count_point, other_configuration[i]);
      }
#pragma omp critical
      for(decltype(size(other_configuration)) i = 0; i < size(other_configuration); ++i) {
        if(dispersion_private[i] != DispersionType{}) {
          dispersion[index_configuration + i] = dispersion_private[i];
        }
      }
}
    } else {
      for(size_t i(0); i < configuration_size; ++i) {
        count_vector[i] = CountType(number_types);
        for(size_t j(0); j < i; ++j) {
          // TODO: varphi(configuration[i], configuration[j]) only needs to be computed once
          AbstractDispersion::template update_count<1>(varphi, count_vector[i][get_type(configuration[j])], configuration[i], configuration[j]);
          AbstractDispersion::template update_count<1>(varphi, count_vector[j][get_type(configuration[i])], configuration[i], configuration[j]);
        }
      }

      if(compute_on_configuration) {
        for(size_t i(0); i < configuration_size; ++i) {
          dispersion[i] = DispersionType(number_types);
          for(size_t j(0); j < i; ++j) {
            // TODO: varphi(configuration[i], configuration[j]) only needs to be computed once
            dispersion[i][get_type(configuration[j])] += AbstractDispersion::template delta<1, true>(varphi, count_vector[j][get_type(configuration[i])], configuration[j], configuration[i]);
            dispersion[j][get_type(configuration[i])] += AbstractDispersion::template delta<1, true>(varphi, count_vector[i][get_type(configuration[j])], configuration[j], configuration[i]);
          }
          add_count_to_dispersion<1, AbstractDispersion, 1>(varphi, dispersion[i], count_vector[i], configuration[i]);
        }
      }

#pragma omp parallel default(none) shared(other_configuration) \
      shared(dispersion, configuration, number_types, count_vector, varphi)
{
      decltype(dispersion) dispersion_private(size(other_configuration));
#pragma omp for nowait
      for(decltype(size(other_configuration)) i = 0; i < size(other_configuration); ++i) {
        dispersion_private[i] = DispersionType(number_types);
        CountType count_point(number_types);
        for(size_t j(0); j < size(configuration); ++j) {
          dispersion_private[i][get_type(configuration[j])] += AbstractDispersion::template delta<1, false>(varphi, count_vector[j][get_type(other_configuration[i])], configuration[j], other_configuration[i]);
          AbstractDispersion::template update_count<0>(varphi, count_point[get_type(configuration[j])], other_configuration[i], configuration[j]);
        }
        add_count_to_dispersion<0, AbstractDispersion, 1>(varphi, dispersion_private[i], count_point, other_configuration[i]);
      }
#pragma omp critical
      for(decltype(size(other_configuration)) i = 0; i < size(other_configuration); ++i) {
        if(dispersion_private[i] != DispersionType{}) {
          dispersion[index_configuration + i] = dispersion_private[i];
        }
      }
}
    }
    return dispersion;
  }
};

} // namespace detail

template<bool ComputeOnConfiguration = true, typename FloatType, typename... Configurations>
inline auto compute_dispersion_for_fitting(const Saturated_model<FloatType>& model,
                                           R_xlen_t number_types,
                                           int nthreads,
                                           Configurations&&... configurations) {
  // TODO: Figure out how to add ComputationType template param to this function
  // TODO: Try to make ComputeOnConfiguration a template parameter, didn't figure out how w/ current version of dispatch
  return detail::dispatch_model<detail::dispersion_computation_fitting>(model, number_types, nthreads, ComputeOnConfiguration, std::forward<Configurations>(configurations)...);
}

} // namespace ppjsdm

#endif // INCLUDE_COMPUTE_DISPERSION_FITTING
