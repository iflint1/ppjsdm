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

// Get the 'Saturation-th' element in the underlying Min-heap with maximum size Saturation + Buffer,
// satisfying a condition.
template<long long int Buffer, long long int Depth = Buffer>
struct get_smallest_if_implementation {
  template<typename Heap, typename Compare, typename Condition>
  auto operator()(unsigned long long int saturation, const Heap& heap, Compare comp, Condition cond) const {
    if(heap.size() == saturation + Buffer - Depth) {
      const auto nth(get_nth<Buffer - Depth>(heap, comp));
      if(!cond(nth)) {
        return nth;
      } else {
        return get_nth<Buffer - Depth + 1>(heap, comp);
      }
    } else {
      return get_smallest_if_implementation<Buffer, Depth - 1>{}(saturation, heap, comp, cond);
    }
  }
};

template<long long int Buffer>
struct get_smallest_if_implementation<Buffer, 0> {
  template<typename Heap, typename Compare, typename Condition>
  auto operator()(unsigned long long int, const Heap& heap, Compare comp, Condition cond) const {
    const auto nth(get_nth<Buffer>(heap, comp));
    if(!cond(nth)) {
      return nth;
    } else {
      return get_nth<Buffer + 1>(heap, comp);
    }
  }
};

template<long long int Buffer, long long int Depth = Buffer, bool IsFirstStep = (Buffer == Depth), bool IsBufferNonZero = (Buffer != 0)>
struct get_nth_smallest_if_implementation;

// Generic iteration of the algorithm
template<long long int Buffer, long long int Depth>
struct get_nth_smallest_if_implementation<Buffer, Depth, false, true> {
  template<typename Heap, typename Compare, typename Condition>
  auto operator()(typename Heap::size_type N, const Heap& heap, Compare comp, Condition cond) const {
    if(heap.size() == N + static_cast<decltype(N)>(Buffer - Depth)) { // In this case, avoid code duplication by running to the last step of the algorithm.
      return get_nth_smallest_if_implementation<Buffer - Depth, 0, false, true>{}(N, heap, comp, cond);
    } else { // Continue iteration
      return get_nth_smallest_if_implementation<Buffer, Depth - 1, false, true>{}(N, heap, comp, cond);
    }
  }
};

// Final step of the generic algorithm
template<long long int Buffer>
struct get_nth_smallest_if_implementation<Buffer, 0, false, true> {
  template<typename Heap, typename Compare, typename Condition>
  auto operator()(typename Heap::size_type, const Heap& heap, Compare comp, Condition cond) const {
    const auto nth(get_nth<Buffer>(heap, comp));
    if(!cond(nth)) {
      return nth;
    } else {
      return get_nth<Buffer - 1>(heap, comp);
    }
  }
};

// When Buffer == 0, Depth is necessarily 0, and we know that (i) we need to return the first element,
// and (ii) it does not satisfy the condition
template<long long int Buffer, long long int Depth>
struct get_nth_smallest_if_implementation<Buffer, Depth, true, false> {
  template<typename Heap, typename Compare, typename Condition>
  auto operator()(typename Heap::size_type, const Heap& heap, Compare comp, Condition) const {
    return get_nth<0>(heap, comp);
  }
};

// When Buffer == Depth and Buffer != 0 we're in the first step of the algorithm,
// and so we want to replicate the generic step of the algorithm,
// noting however that the condition is known to not be satisfied.
template<long long int Buffer, long long int Depth>
struct get_nth_smallest_if_implementation<Buffer, Depth, true, true> {
  template<typename Heap, typename Compare, typename Condition>
  auto operator()(typename Heap::size_type N, const Heap& heap, Compare comp, Condition cond) const {
    if(heap.size() == N) {
      return get_nth<0>(heap, comp);
    } else { // Continue iteration
      return get_nth_smallest_if_implementation<Buffer, Depth - 1, false, true>{}(N, heap, comp, cond);
    }
  }
};

} // namespace detail

// Given a heap x_1 >  ... > x_k, and k >= N, returns x_{1 + k - N} if !cond(x_{1 + k - N}), and if not x_{k - N}.
// In words, this returns the N-th smallest element in the heap if it does not satisfy a condition, and the 'N+1'-th otherwise.
// It is known at compile-time that 0 <= k - N <= Buffer.
// It is additionally known at compile-time that if k == N, then the element does not satisfy the condition.
template<long long int Buffer, typename Heap, typename Compare, typename Condition>
inline auto get_nth_smallest_if(typename Heap::size_type N, const Heap& heap, Compare comp, Condition cond) {
  return detail::get_nth_smallest_if_implementation<Buffer>{}(N, heap, comp, cond);
}

template<long long int Buffer, typename Heap, typename Condition>
inline auto get_nth_smallest_if(typename Heap::size_type N, const Heap& heap, Condition cond) {
  return get_nth_smallest_if<Buffer>(N, heap, std::less<std::remove_reference_t<std::remove_cv_t<decltype(heap[0])>>>{}, cond);
}

template<long long int Buffer, typename Vector, typename Compare>
inline auto get_smallest(unsigned long long int saturation, const Vector& vector, Compare comp) {
  return detail::get_smallest_if_implementation<Buffer>{}(saturation, vector, comp, [](const auto){ return false; });
}

template<long long int N, typename Vector>
inline auto get_smallest(unsigned long long int saturation, const Vector& vector) {
  return get_smallest<N>(saturation, vector, std::less<std::remove_reference_t<std::remove_cv_t<decltype(vector[0])>>>{});
}

template<long long int Buffer, typename Vector, typename Compare>
inline auto get_smallest_excluding(unsigned long long int saturation,
                                   const Vector& vector,
                                   const typename Vector::value_type& value,
                                   Compare comp) {
  return detail::get_smallest_if_implementation<Buffer>{}(saturation, vector, comp, [value](const auto& element){ return element == value; });
}

template<long long int N, typename Vector>
inline auto get_smallest_excluding(unsigned long long int saturation,
                                   const Vector& vector,
                                   const typename Vector::value_type& value) {
  return get_smallest_excluding<N>(saturation,
                                   vector,
                                   value,
                                   std::less<std::remove_reference_t<std::remove_cv_t<decltype(vector[0])>>>{});
}

} // namespace ppjsdm

#endif // INCLUDE_HEAP
