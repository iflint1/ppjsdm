#ifndef INCLUDE_PPJSDM_FOR_EACH_CONTAINER
#define INCLUDE_PPJSDM_FOR_EACH_CONTAINER

#include <utility> // std::forward

namespace ppjsdm {

template<typename... Args>
void for_each_container(Args&&...) {
  return;
}

// Apply a function over a variadic number of containers which can be indexed but do not necessarily have iterators.
template<typename F, typename Container, typename... Others>
void for_each_container(F&& f, const Container& container, Others&&... others) {
  for(const auto& element: container) {
    f(element);
  }
  for_each_container(std::forward<F>(f), std::forward<Others>(others)...);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_FOR_EACH_CONTAINER
