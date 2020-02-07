#ifndef INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST
#define INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST

#include <Rinternals.h>

#include "../utility/backwards_markov_chain.hpp"

namespace ppjsdm {

template<typename Configuration, typename Model>
inline auto simulate_coupling_from_the_past(const Model& model, R_xlen_t number_types) {
  auto Z(make_backwards_markov_chain<Configuration>(model, number_types));
  Z.extend_until_T0();
  while(true) {
    const auto coalescence(Z.compute_LU_and_check_coalescence());
    if(coalescence.first) {
      return coalescence.second;
    }
    R_CheckUserInterrupt();
    Z.extend_backwards(Z.size());
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST
