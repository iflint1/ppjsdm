#ifndef INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST
#define INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST

#include <RcppThread.h>
#include <Rinternals.h>

#include "../utility/backwards_markov_chain.hpp"

namespace ppjsdm {

template<typename Configuration, typename Generator, typename Model>
inline auto simulate_coupling_from_the_past(Generator& generator,
                                            const Model& model) {
  auto Z(make_backwards_markov_chain<Configuration>(model, generator));
  Z.extend_until_T0(generator);
  while(true) {
    const auto coalescence(Z.compute_LU_and_check_coalescence(generator));
    if(coalescence.first) {
      return coalescence.second;
    }
    RcppThread::checkUserInterrupt();
    Z.extend_backwards(generator, Z.size());
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_COUPLING_FROM_THE_PAST
