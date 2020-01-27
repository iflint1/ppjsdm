#ifndef INCLUDE_PPJSDM_SIZE_T
#define INCLUDE_PPJSDM_SIZE_T

#include <utility> // std::declval

namespace ppjsdm {

template<typename T>
using size_t = decltype(size(std::declval<T>()));

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SIZE_T
