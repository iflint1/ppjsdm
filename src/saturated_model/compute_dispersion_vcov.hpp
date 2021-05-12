#ifndef INCLUDE_COMPUTE_DISPERSION_VCOV
#define INCLUDE_COMPUTE_DISPERSION_VCOV

#include <Rcpp.h>

#include "compute_dispersion_implementation.hpp"
#include "saturated_model.hpp"
#include "../utility/flatten_strict_upper_triangular.hpp"

#include <tuple> // std::pair
#include <type_traits> // std::remove_cv_t
#include <utility> // std::forward
#include <vector> // std::vector

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ppjsdm {
namespace detail {

template<typename AbstractDispersion>
struct generic_vcov_dispersion_computation {
  template<typename Configuration>
  auto operator()(const Saturated_model& varphi,
                R_xlen_t number_types,
                const Configuration& configuration,
                const Configuration& restricted_configuration,
                typename Configuration::size_type min_index,
                typename Configuration::size_type max_index,
                int nthreads) const {
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif

    // TODO: Better way?
    min_index = std::max<typename Configuration::size_type>(min_index, 0);
    max_index = std::min<typename Configuration::size_type>(max_index, ppjsdm::size(restricted_configuration) * (ppjsdm::size(restricted_configuration) - 1) / 2);

    using ValueType = typename AbstractDispersion::ValueType;
    using CountType = std::vector<ValueType>;
    using DispersionType = std::vector<double>;

    const auto restricted_configuration_size(size(restricted_configuration));
    const auto configuration_size(size(configuration));
    using size_t = std::remove_cv_t<decltype(size(configuration))>;
    std::vector<CountType> count_vector(configuration_size);

    std::vector<size_t> index_in_configuration(restricted_configuration_size);
    for(size_t i(0); i < restricted_configuration_size; ++i) {
      for(size_t j(0); j < configuration_size; ++j) {
        if(is_equal(restricted_configuration[i], configuration[j])) {
          index_in_configuration[i] = j;
          break;
        }
      }
    }

    // TODO: In fact, we can also run this case when saturation >= max(size_by_type)
    if(static_cast<decltype(size(configuration))>(varphi.get_saturation()) >= configuration_size) {
      for(size_t i(0); i < configuration_size; ++i) {
        count_vector[i] = CountType(number_types);
        for(size_t j(0); j < i; ++j) {
          // TODO: varphi(configuration[i], configuration[j]) only needs to be computed once
          AbstractDispersion::template update_count<std::numeric_limits<int>::max()>(varphi, count_vector[i][get_type(configuration[j])], configuration[i], configuration[j]);
          AbstractDispersion::template update_count<std::numeric_limits<int>::max()>(varphi, count_vector[j][get_type(configuration[i])], configuration[i], configuration[j]);
        }
      }

      std::vector<DispersionType> dispersion_i(max_index - min_index);
      std::vector<DispersionType> dispersion_j(dispersion_i.size());

      for(std::remove_cv_t<decltype(dispersion_i.size())> index(min_index); index < max_index; ++index) {
        const auto pr(decode_linear(index, restricted_configuration_size));
        const size_t i(pr.first);
        const size_t j(pr.second);

        dispersion_i[index - min_index] = DispersionType(number_types);
        dispersion_j[index - min_index] = DispersionType(number_types);

        add_count_to_dispersion_discarding<0, AbstractDispersion, 2>(varphi, dispersion_i[index - min_index],
                                                                     count_vector[index_in_configuration[i]], restricted_configuration[i], restricted_configuration[j]);
        add_count_to_dispersion_discarding<0, AbstractDispersion, 2>(varphi, dispersion_j[index - min_index],
                                                                     count_vector[index_in_configuration[j]], restricted_configuration[j], restricted_configuration[i]);
      }

      return std::pair<std::vector<DispersionType>, std::vector<DispersionType>>(dispersion_i, dispersion_j);
    } else {
      for(size_t i(0); i < configuration_size; ++i) {
        count_vector[i] = CountType(number_types);
        for(size_t j(0); j < i; ++j) {
          // TODO: varphi(configuration[i], configuration[j]) only needs to be computed once
          AbstractDispersion::template update_count<2>(varphi, count_vector[i][get_type(configuration[j])], configuration[i], configuration[j]);
          AbstractDispersion::template update_count<2>(varphi, count_vector[j][get_type(configuration[i])], configuration[i], configuration[j]);
        }
      }

      std::vector<DispersionType> sum_deltas;
      if(number_types > 1) {
        sum_deltas = std::vector<DispersionType>(configuration_size);
        for(size_t i(0); i < configuration_size; ++i) {
          sum_deltas[i] = DispersionType(number_types);
          for(size_t j(0); j < i; ++j) {
            sum_deltas[i][get_type(configuration[j])] += AbstractDispersion::template delta<2, true>(varphi, count_vector[j][get_type(configuration[i])], configuration[j], configuration[i]);
            sum_deltas[j][get_type(configuration[i])] += AbstractDispersion::template delta<2, true>(varphi, count_vector[i][get_type(configuration[j])], configuration[i], configuration[j]);
          }
        }
      }


      std::vector<DispersionType> dispersion_i(max_index - min_index);
      std::vector<DispersionType> dispersion_j(dispersion_i.size());

#pragma omp parallel default(none) shared(restricted_configuration) \
      shared(min_index, max_index, count_vector, dispersion_i, dispersion_j) \
      shared(number_types, index_in_configuration, varphi, sum_deltas)
{
      decltype(dispersion_i) dispersion_i_private(dispersion_i.size());
      decltype(dispersion_j) dispersion_j_private(dispersion_i.size());
#pragma omp for nowait
      for(std::remove_cv_t<decltype(dispersion_i.size())> index = min_index; index < max_index; ++index) {
        const auto pr(decode_linear(index, restricted_configuration_size));
        const size_t i(pr.first);
        const size_t j(pr.second);

        // At this point, we have i < j < restricted_configuration_size
        if(get_type(restricted_configuration[i]) == get_type(restricted_configuration[j])) {
          dispersion_i_private[index - min_index] = DispersionType(number_types);
          dispersion_j_private[index - min_index] = DispersionType(number_types);
          for(size_t l(0); l < restricted_configuration_size; ++l) {
            if(l != i && l != j) {
              dispersion_i_private[index - min_index][get_type(restricted_configuration[l])] += AbstractDispersion::template delta_discarding<2, true>(varphi, count_vector[index_in_configuration[l]][get_type(restricted_configuration[i])], restricted_configuration[l], restricted_configuration[i], restricted_configuration[j]);
              dispersion_j_private[index - min_index][get_type(restricted_configuration[l])] += AbstractDispersion::template delta_discarding<2, true>(varphi, count_vector[index_in_configuration[l]][get_type(restricted_configuration[j])], restricted_configuration[l], restricted_configuration[j], restricted_configuration[i]);
            }
          }
        } else {
          dispersion_i_private[index - min_index] = sum_deltas[index_in_configuration[i]];
          dispersion_j_private[index - min_index] = sum_deltas[index_in_configuration[j]];

          dispersion_i_private[index - min_index][get_type(restricted_configuration[j])] -= AbstractDispersion::template delta<2, true>(varphi, count_vector[index_in_configuration[j]][get_type(restricted_configuration[i])], restricted_configuration[j], restricted_configuration[i]);
          dispersion_j_private[index - min_index][get_type(restricted_configuration[i])] -= AbstractDispersion::template delta<2, true>(varphi, count_vector[index_in_configuration[i]][get_type(restricted_configuration[j])], restricted_configuration[i], restricted_configuration[j]);
        }
        add_count_to_dispersion_discarding<2, AbstractDispersion, 1>(varphi, dispersion_i_private[index - min_index],
                                                                     count_vector[index_in_configuration[i]], restricted_configuration[i], restricted_configuration[j]);
        add_count_to_dispersion_discarding<2, AbstractDispersion, 1>(varphi, dispersion_j_private[index - min_index],
                                                                     count_vector[index_in_configuration[j]], restricted_configuration[j], restricted_configuration[i]);
      }
#pragma omp critical
      for(std::remove_cv_t<decltype(dispersion_i.size())> index(min_index); index < max_index; ++index) {
        if(dispersion_i_private[index - min_index] != DispersionType{}) {
          dispersion_i[index - min_index] = dispersion_i_private[index - min_index];
          dispersion_j[index - min_index] = dispersion_j_private[index - min_index];
        }
      }
}
      return std::pair<std::vector<DispersionType>, std::vector<DispersionType>>(dispersion_i, dispersion_j);
    }
  }
};

} // namespace detail

template<typename Configuration>
inline auto compute_dispersion_for_vcov(const Saturated_model& model,
                                        R_xlen_t number_types,
                                        const Configuration& configuration,
                                        const Configuration& restricted_configuration,
                                        typename Configuration::size_type min_index,
                                        typename Configuration::size_type max_index,
                                        int nthreads = 1) {
  return detail::dispatch_model<detail::generic_vcov_dispersion_computation>(model, number_types, configuration, restricted_configuration, min_index, max_index, nthreads);
}

template<typename Configuration>
inline auto compute_dispersion_for_vcov(const Saturated_model& model,
                                        R_xlen_t number_types,
                                        const Configuration& configuration,
                                        int nthreads = 1) {
  return compute_dispersion_for_vcov(model, number_types, configuration, configuration, 0, size(configuration) * (size(configuration) - 1) / 2, nthreads);
}

} // namespace ppjsdm

#endif // INCLUDE_COMPUTE_DISPERSION_VCOV
