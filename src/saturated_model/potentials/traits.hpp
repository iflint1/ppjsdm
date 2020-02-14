#ifndef INCLUDE_PPJSDM_POTENTIALS_TRAITS
#define INCLUDE_PPJSDM_POTENTIALS_TRAITS

#include <type_traits> // std::enable_if, std::false_type, std::true_type
#include <utility> // std::declval

namespace ppjsdm {
// TODO: Write a couple of static_asserts for type traits below.

template<typename T, typename = void>
struct has_nonzero_value: std::false_type {};

template<typename T>
struct has_nonzero_value<T, decltype(static_cast<void>(T::nonzero_value), void())>: std::true_type {};

template<typename T>
constexpr bool has_nonzero_value_v = has_nonzero_value<T>::value;

template<typename T, typename = void>
struct has_square_lower_endpoint: std::false_type {};

template<typename T>
struct has_square_lower_endpoint<T, decltype(static_cast<void>(std::declval<T>().get_square_lower_endpoint(0.)), void())>: std::true_type {};

template<typename T>
constexpr bool has_square_lower_endpoint_v = has_square_lower_endpoint<T>::value;

template<typename T, typename = void>
struct is_nonincreasing_after_lower_endpoint: std::false_type {};

template<typename T>
struct is_nonincreasing_after_lower_endpoint<T, std::enable_if_t<has_square_lower_endpoint_v<T>>> {
  static constexpr bool value = T::is_nonincreasing_after_lower_endpoint;
};

template<typename T>
constexpr bool is_nonincreasing_after_lower_endpoint_v = is_nonincreasing_after_lower_endpoint<T>::value;

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_POTENTIALS_TRAITS
