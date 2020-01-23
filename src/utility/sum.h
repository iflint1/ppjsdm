#ifndef INCLUDE_PPJSDM_SUM
#define INCLUDE_PPJSDM_SUM

#include "../utility/size_t.h"

namespace ppjsdm {

// Note: Can't use std::accumulate instead of this because of ``ambiguous overload'' when compiling on Windows.
// Note: cannot deduce the return type since I want this to also work on lists whose underlying type is not known at compile time.
template<typename ReturnType, typename Vector>
inline auto sum(const Vector& vector, size_t<Vector> size) {
  ReturnType sum(0);
  for(decltype(size) i(0); i < size; ++i) {
    sum += vector[i];
  }
  return sum;
}

template<typename ReturnType, typename Vector>
inline auto sum(const Vector& vector) {
  return sum<ReturnType>(vector, vector.size());
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SUM
