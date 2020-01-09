#ifndef INCLUDE_PPJSDM_SUM
#define INCLUDE_PPJSDM_SUM

#include <Rcpp.h>

namespace ppjsdm {

// Note: cannot deduce the return type since I want this to also work on lists whose type is not known at compile time.
template<typename ReturnType, typename T>
inline auto sum(const T& vector, R_xlen_t size) {
  ReturnType sum(0);
  for(R_xlen_t i(0); i < size; ++i) {
    sum += vector[i];
  }
  return sum;
}

template<typename ReturnType, typename T>
inline auto sum(const T& vector) {
  return sum<ReturnType>(vector, vector.size());
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SUM
