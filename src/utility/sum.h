#ifndef INCLUDE_PPJSDM_SUM
#define INCLUDE_PPJSDM_SUM

#include <utility> // std::declval

namespace ppjsdm {
namespace detail {

template<typename T>
using size_t = decltype(std::declval<T>().size());

} // namespace detail

// Note: cannot deduce the return type since I want this to also work on lists whose underlying type is not known at compile time.
template<typename ReturnType, typename T>
inline auto sum(const T& vector_or_list, detail::size_t<T> size) {
  ReturnType sum(0);
  for(decltype(size) i(0); i < size; ++i) {
    sum += vector_or_list[i];
  }
  return sum;
}

template<typename ReturnType, typename T>
inline auto sum(const T& vector_or_list) {
  return sum<ReturnType>(vector_or_list, vector_or_list.size());
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SUM
