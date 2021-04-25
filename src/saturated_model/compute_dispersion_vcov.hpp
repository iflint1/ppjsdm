#ifndef INCLUDE_COMPUTE_DISPERSION
#define INCLUDE_COMPUTE_DISPERSION

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
                int nthreads) const {
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif

    using ValueType = typename AbstractDispersion::ValueType;
    using CountType = std::vector<ValueType>;
    using DispersionType = std::vector<double>;

    const auto configuration_size(size(configuration));
    using size_t = std::remove_cv_t<decltype(size(configuration))>;
    std::vector<CountType> count_vector(configuration_size);

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

      std::vector<DispersionType> dispersion_i(configuration_size * (configuration_size - 1) / 2);
      std::vector<DispersionType> dispersion_j(dispersion_i.size());

      for(size_t i(0); i < configuration_size; ++i) {
        for(size_t j(i + 1); j < configuration_size; ++j) {
          const auto index(encode_linear(i, j, configuration_size));

          dispersion_i[index] = DispersionType(number_types);
          dispersion_j[index] = DispersionType(number_types);

          add_count_to_dispersion_discarding<0, AbstractDispersion, 2>(varphi, dispersion_i[index],
                                                                       count_vector[i], configuration[i], configuration[j]);
          add_count_to_dispersion_discarding<0, AbstractDispersion, 2>(varphi, dispersion_j[index],
                                                                       count_vector[j], configuration[j], configuration[i]);
        }
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


      std::vector<DispersionType> dispersion_i(configuration_size * (configuration_size - 1) / 2);
      std::vector<DispersionType> dispersion_j(dispersion_i.size());

#pragma omp parallel for
      for(std::remove_cv_t<decltype(dispersion_i.size())> index = 0; index < dispersion_i.size(); ++index) {
        const auto pr(decode_linear(index, configuration_size));
        const size_t i(pr.first);
        const size_t j(pr.second);

        // At this point, we have i < j < configuration_size
        if(get_type(configuration[i]) == get_type(configuration[j])) {
          dispersion_i[index] = DispersionType(number_types);
          dispersion_j[index] = DispersionType(number_types);
          for(size_t l(0); l < configuration_size; ++l) {
            if(l != i && l != j) {
              dispersion_i[index][get_type(configuration[l])] += AbstractDispersion::template delta_discarding<2, true>(varphi, count_vector[l][get_type(configuration[i])], configuration[l], configuration[i], configuration[j]);
              dispersion_j[index][get_type(configuration[l])] += AbstractDispersion::template delta_discarding<2, true>(varphi, count_vector[l][get_type(configuration[j])], configuration[l], configuration[j], configuration[i]);
            }
          }
        } else {
          dispersion_i[index] = sum_deltas[i];
          dispersion_j[index] = sum_deltas[j];

          dispersion_i[index][get_type(configuration[j])] -= AbstractDispersion::template delta<2, true>(varphi, count_vector[j][get_type(configuration[i])], configuration[j], configuration[i]);
          dispersion_j[index][get_type(configuration[i])] -= AbstractDispersion::template delta<2, true>(varphi, count_vector[i][get_type(configuration[j])], configuration[i], configuration[j]);
        }
        add_count_to_dispersion_discarding<2, AbstractDispersion, 1>(varphi, dispersion_i[index],
                                                                     count_vector[i], configuration[i], configuration[j]);
        add_count_to_dispersion_discarding<2, AbstractDispersion, 1>(varphi, dispersion_j[index],
                                                                     count_vector[j], configuration[j], configuration[i]);
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
                                        int nthreads) {
  return detail::dispatch_model<detail::generic_vcov_dispersion_computation>(model, number_types, configuration, nthreads);
}

} // namespace ppjsdm

#endif // INCLUDE_COMPUTE_DISPERSION
