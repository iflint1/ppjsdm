#ifndef INCLUDE_FIND_IN_CONTAINERS
#define INCLUDE_FIND_IN_CONTAINERS

#include <algorithm> // std::find_if
#include <utility> // std::forward

namespace ppjsdm {

template<typename UnaryPredicate>
inline auto find_in_containers(UnaryPredicate&&) {
  return false;
}

// Is a point in the variadic containers?
template<typename UnaryPredicate, typename Container, typename... Others>
inline auto find_in_containers(UnaryPredicate&& p, const Container& container, Others&&... others) {
  return (std::find_if(container.begin(), container.end(), p) != container.end()) || find_in_containers(std::forward<UnaryPredicate>(p), std::forward<Others>(others)...);
}

} // namespace ppjsdm

#endif // INCLUDE_FIND_IN_CONTAINERS
