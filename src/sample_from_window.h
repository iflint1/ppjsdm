#ifndef INCLUDE_PPJSDM_SAMPLE_FROM_WINDOW
#define INCLUDE_PPJSDM_SAMPLE_FROM_WINDOW

#include <Rinternals.h>

namespace ppjsdm {

template<typename S, typename T, typename U, typename V>
inline void sample_from_window(const S& window, T& x, U& y, V& types, R_xlen_t additions, R_xlen_t filling, int types_value) {
  for(R_xlen_t i{0}; i < additions; ++i) {
    const auto sample{window.sample()};
    x[filling + i] = sample.first;
    y[filling + i] = sample.second;
    types[filling + i] = types_value;
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SAMPLE_FROM_WINDOW
