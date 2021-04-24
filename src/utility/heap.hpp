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
    if(count.size() > 2 && comp(count[1], count[2])) {
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
    if(comp(count[1], count[2])) {
      if(count.size() > 6) {
        if(comp(count[1], count[5])) {
          if(comp(count[5], count[6])) {
            return count[6];
          } else {
            return count[5];
          }
        } else {
          if(comp(count[1], count[6])) {
            return count[6];
          } else {
            return count[1];
          }
        }
      } else if(count.size() == 6 && comp(count[1], count[5])) {
        return count[5];
      } else {
        return count[1];
      }
    } else {
      if(count.size() > 4) {
        if(comp(count[2], count[3])) {
          if(comp(count[3], count[4])) {
            return count[4];
          } else {
            return count[3];
          }
        } else {
          if(comp(count[2], count[4])) {
            return count[4];
          } else {
            return count[2];
          }
        }
      } else if(count.size() == 4 && comp(count[2], count[3])) {
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

namespace detail {
// Get the 'Saturation-th' element in the underlying Min-heap with maximum size Saturation + Buffer.
template<long long int Buffer, long long int Depth = Buffer>
struct get_smallest_implementation {
  template<typename Vector, typename Compare>
  auto operator()(unsigned long long int saturation, const Vector& count, Compare comp) const {
    if(count.size() == saturation + Buffer - Depth) {
      // Note that since saturation >= 1, the calls to get_nth satisfy count.size() >= N + 1
      return get_nth<Buffer - Depth>(count, comp);
    } else {
      return get_smallest_implementation<Buffer, Depth - 1>{}(saturation, count, comp);
    }
  }
};

template<long long int Buffer>
struct get_smallest_implementation<Buffer, 0> {
  template<typename Vector, typename Compare>
  auto operator()(unsigned long long int, const Vector& count, Compare comp) const {
    return get_nth<Buffer>(count, comp);
  }
};

// Get the 'Saturation-th' element in the underlying Min-heap with maximum size Saturation + Buffer,
// excluding an element.
template<long long int Buffer, long long int Depth = Buffer>
struct get_smallest_excluding_implementation {
  template<typename Vector, typename Compare>
  auto operator()(unsigned long long int saturation, const Vector& count, const typename Vector::value_type& value, Compare comp) const {
    if(count.size() == saturation + Buffer - Depth) {
      const auto nth(get_nth<Buffer - Depth>(count, comp));
      if(nth != value) {
        return nth;
      } else {
        return get_nth<Buffer - Depth + 1>(count, comp);
      }
    } else {
      return get_smallest_excluding_implementation<Buffer, Depth - 1>{}(saturation, count, value, comp);
    }
  }
};

template<long long int Buffer>
struct get_smallest_excluding_implementation<Buffer, 0> {
  template<typename Vector, typename Compare>
  auto operator()(unsigned long long int, const Vector& count, const typename Vector::value_type& value, Compare comp) const {
    const auto nth(get_nth<Buffer>(count, comp));
    if(nth != value) {
      return nth;
    } else {
      return get_nth<Buffer + 1>(count, comp);
    }
  }
};

} // namespace detail

template<long long int Buffer, typename Vector, typename Compare>
inline auto get_smallest(unsigned long long int saturation, const Vector& vector, Compare comp) {
  return detail::get_smallest_implementation<Buffer>{}(saturation, vector, comp);
}

template<long long int N, typename Vector>
inline auto get_smallest(unsigned long long int saturation, const Vector& vector) {
  return get_smallest<N>(saturation, vector, std::less<std::remove_reference_t<std::remove_cv_t<decltype(vector[0])>>>{});
}

template<long long int Buffer, typename Vector, typename Compare>
inline auto get_smallest_excluding(unsigned long long int saturation, const Vector& vector, const typename Vector::value_type& value, Compare comp) {
  return detail::get_smallest_excluding_implementation<Buffer>{}(saturation, vector, value, comp);
}

template<long long int N, typename Vector>
inline auto get_smallest_excluding(unsigned long long int saturation, const Vector& vector, const typename Vector::value_type& value) {
  return get_smallest_excluding<N>(saturation, vector, value, std::less<std::remove_reference_t<std::remove_cv_t<decltype(vector[0])>>>{});
}

} // namespace ppjsdm

#endif // INCLUDE_HEAP
