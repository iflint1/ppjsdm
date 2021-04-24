#ifndef INCLUDE_HEAP
#define INCLUDE_HEAP

#include <type_traits> // std::remove_reference_t, std::remove_cv_t

namespace ppjsdm {
namespace detail {

// Note: In general, this function assumes that count.size() >= N + 1
template<long long int N>
struct get_nth_implementation {
  template<typename Vector, typename Compare>
  auto operator()(const Vector& count, Compare comp) const;
};

// Note: Assumes that count.size() >= 1
template<>
struct get_nth_implementation<0> {
  template<typename Vector, typename Compare>
  auto operator()(const Vector& count, Compare) const {
    return count[0];
  }
};

// Note: Assumes that count.size() >= 2
template<>
struct get_nth_implementation<1> {
  template<typename Vector, typename Compare>
  auto operator()(const Vector& count, Compare comp) const {
    if(count.size() > 2 && comp(count[2], count[1])) {
      return count[2];
    } else {
      return count[1];
    }
  }
};

// Note: Assumes that count.size() >= 3
template<>
struct get_nth_implementation<2> {
  template<typename Vector, typename Compare>
  auto operator()(const Vector& count, Compare comp) const {
    if(comp(count[2], count[1])) {
      if(count.size() > 6) {
        if(comp(count[5], count[1])) {
          if(comp(count[6], count[5])) {
            return count[6];
          } else {
            return count[5];
          }
        } else {
          if(comp(count[6], count[1])) {
            return count[6];
          } else {
            return count[1];
          }
        }
      } else if(count.size() == 6 && comp(count[5], count[1])) {
        return count[5];
      } else {
        return count[1];
      }
    } else {
      if(count.size() > 4) {
        if(comp(count[3], count[2])) {
          if(comp(count[4], count[3])) {
            return count[4];
          } else {
            return count[3];
          }
        } else {
          if(comp(count[4], count[2])) {
            return count[4];
          } else {
            return count[2];
          }
        }
      } else if(count.size() == 4 && comp(count[3], count[2])) {
        return count[3];
      } else {
        return count[2];
      }
    }
  }
};

} // namespace detail

template<long long int N, typename Vector, typename Compare>
inline auto get_nth(const Vector& vector, Compare comp) {
  return detail::get_nth_implementation<N>{}(vector, comp);
}

template<long long int N, typename Vector>
inline auto get_nth(const Vector& vector) {
  return get_nth<N>(vector, std::less<std::remove_reference_t<std::remove_cv_t<decltype(vector[0])>>>{});
}

} // namespace ppjsdm

#endif // INCLUDE_HEAP
