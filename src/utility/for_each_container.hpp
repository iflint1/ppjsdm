#ifndef INCLUDE_PPJSDM_FOR_EACH_CONTAINER
#define INCLUDE_PPJSDM_FOR_EACH_CONTAINER

#include <utility> // std::forward

namespace ppjsdm {

template<typename F, typename Condition>
inline void conditional_for_each_container(F&&, Condition&&) {
  return;
}

// Apply a function over a variadic number of containers which can be indexed but do not necessarily have iterators.
// Stop the loop when condition is satisfied.
template<typename F, typename Condition, typename Container, typename... Others>
inline void conditional_for_each_container(F&& f, Condition&& condition, const Container& container, Others&&... others) {
  for(const auto& element: container) {
    f(element);
    if(condition()) {
      return;
    }
  }
  conditional_for_each_container(std::forward<F>(f), std::forward<Condition>(condition), std::forward<Others>(others)...);
}

// Apply a function over a variadic number of containers which can be indexed but do not necessarily have iterators.
template<typename F, typename... Containers>
inline void for_each_container(F&& f,  Containers&&... containers) {
  conditional_for_each_container(std::forward<F>(f), [](){ return false; }, std::forward<Containers>(containers)...);
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_FOR_EACH_CONTAINER
