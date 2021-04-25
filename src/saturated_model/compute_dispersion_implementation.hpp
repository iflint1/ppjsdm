#ifndef INCLUDE_COMPUTE_DISPERSION_IMPLEMENTATION
#define INCLUDE_COMPUTE_DISPERSION_IMPLEMENTATION

#include <Rcpp.h>

#include "saturated_model.hpp"
#include "../point/point_manipulation.hpp"
#include "../point/square_distance.hpp"
#include "../utility/heap.hpp"

#include <algorithm> // std::accumulate, std::pop_heap, std::push_heap
#include <functional> // std::greater
#include <utility> // std::forward, std::move
#include <vector> // std::vector

namespace ppjsdm {
namespace detail {

template<typename Point, typename Other>
inline auto apply_potential(const Saturated_model& potential, const Point& point, const Other& other) {
  return potential.apply(normalized_square_distance(point, other), get_type(point), get_type(other));
}

struct compute_dispersion_two_values_implementation {
  using ValueType = unsigned long long int;

  template<int Buffer, typename Vector, typename Point>
  static auto add_count_to_dispersion(const Saturated_model& varphi,
                                      const Vector& count_vector,
                                      const Point&,
                                      decltype(count_vector.size()) i) {
    return std::min(varphi.get_saturation(), count_vector[i]);
  }

  template<int Buffer, typename Point, typename ToDiscard, typename Vector>
  static auto add_count_to_dispersion_discarding(const Saturated_model& varphi,
                                                 const Vector& count_vector,
                                                 const Point& point,
                                                 const ToDiscard& to_discard,
                                                 decltype(count_vector.size()) i) {
    return std::min(varphi.get_saturation(), count_vector[i] - static_cast<int>(apply_potential(varphi, point, to_discard) > 0.));
  }

  template<int Buffer = 1, typename Point, typename Other>
  static void update_count(const Saturated_model& varphi, ValueType& count, const Point& point, const Other& other) {
    if((count < varphi.get_saturation() + Buffer) && (apply_potential(varphi, point, other) > 0.)) {
      ++count;
    }
  }

  template<int Buffer, bool CountedAlready, typename Point, typename Other>
  static auto delta(const Saturated_model& varphi,
                    const ValueType& count,
                    const Point& point,
                    const Other& other) {
    if((count < varphi.get_saturation() + static_cast<int>(CountedAlready)) && (apply_potential(varphi, point, other) > 0.)) {
      return 1.;
    } else {
      return 0.;
    }
  }

  template<int Buffer, bool CountedAlready, typename Point, typename Other, typename ToDiscard>
  static auto delta_discarding(const Saturated_model& varphi,
                               const ValueType& count,
                               const Point& point,
                               const Other& other,
                               const ToDiscard&) {
    return delta<Buffer, CountedAlready>(varphi, count - 1, point, other);
  }
};

struct put_potential {
  using HeapOrder = std::greater<double>;

  template<typename Point, typename Other>
  static auto put(const Saturated_model& varphi, double square_distance, const Point& point, const Other& other) {
    return varphi.apply(square_distance, get_type(point), get_type(other));
  }

  static auto get(const Saturated_model&, double element, int, int) {
    return element;
  }
};

struct put_square_distance {
  using HeapOrder = std::less<double>;

  template<typename Point, typename Other>
  static auto put(const Saturated_model&, double square_distance, const Point&, const Other&) {
    return square_distance;
  }

  static auto get(const Saturated_model& varphi, double element, int type, int other_type) {
    return varphi.apply(element, type, other_type);
  }
};

template<typename PutInHeap = put_potential>
class compute_dispersion_generic_implementation {
private:
  using HeapOrder = typename PutInHeap::HeapOrder;
public:
  using ValueType = std::vector<double>;

  template<int Buffer, typename Point, typename Vector>
  static auto add_count_to_dispersion(const Saturated_model& varphi,
                                      const Vector& count_vector,
                                      const Point& point,
                                      decltype(count_vector.size()) i) {
    return accumulate_n_smallest<Buffer>(varphi.get_saturation(), count_vector[i], HeapOrder{}, [&varphi, i, &point](auto count, auto val) {
      return count + PutInHeap::get(varphi, val, i, get_type(point));
    });
  }

  template<int Buffer, typename Point, typename ToDiscard, typename Vector>
  static auto add_count_to_dispersion_discarding(const Saturated_model& varphi,
                                                 const Vector& count_vector,
                                                 const Point& point,
                                                 const ToDiscard& to_discard,
                                                 decltype(count_vector.size()) i) {
    // const auto discard_sq(normalized_square_distance(point, to_discard));
    // const auto discard_disp(PutInHeap::put(varphi, discard_sq, point, to_discard));
    // Vector count_vector_copy(count_vector);
    // count_vector_copy[i].erase(std::remove(count_vector_copy[i].begin(), count_vector_copy[i].end(), discard_disp), count_vector_copy[i].end());
    // std::make_heap(count_vector_copy[i].begin(), count_vector_copy[i].end(), HeapOrder{});
    // return add_count_to_dispersion(varphi, count_vector_copy, point, i);

    const auto discard_sq(normalized_square_distance(point, to_discard));
    if(discard_sq >= varphi.get_square_lower_endpoint(get_type(point), get_type(to_discard))) {
      const auto smallest(get_nth<0>(count_vector[i], HeapOrder{}));
      const auto discard_disp(PutInHeap::put(varphi, discard_sq, point, to_discard));
      if(HeapOrder{}(discard_disp, smallest)) { // to_discard is in count
        return accumulate_n_smallest<Buffer>(varphi.get_saturation() + 1, count_vector[i], HeapOrder{}, [discard_disp, &varphi, i, &point](auto count, auto val) {
          if(val != discard_disp) {
            return count + PutInHeap::get(varphi, val, i, get_type(point));
          } else {
            return count;
          }
        });
      }
    }
    return add_count_to_dispersion<Buffer>(varphi, count_vector, point, i);
  }

  template<int Buffer, typename Point, typename Other>
  static void update_count(const Saturated_model& varphi, ValueType& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(sq >= varphi.get_square_lower_endpoint(get_type(point), get_type(other))) {
      const auto disp(PutInHeap::put(varphi, sq, point, other));
      if(count.size() < varphi.get_saturation() + Buffer) {
        count.emplace_back(disp);
        std::push_heap(count.begin(), count.end(), HeapOrder{});
      } else if(HeapOrder{}(disp, get_nth<0>(count, HeapOrder{}))) {
        count.emplace_back(disp);
        std::pop_heap(count.begin(), count.end(), HeapOrder{});
        count.pop_back();
      }
    }
  }

  template<int Buffer, bool CountedAlready, typename Point, typename Other>
  static auto delta(const Saturated_model& varphi,
                    const ValueType& count,
                    const Point& point,
                    const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(sq >= varphi.get_square_lower_endpoint(get_type(point), get_type(other))) {
      if(count.size() < varphi.get_saturation() + static_cast<int>(CountedAlready)) {
        return varphi.apply(sq, get_type(point), get_type(other));
      } else {
        const auto disp(PutInHeap::put(varphi, sq, point, other));
        const auto smallest(get_nth_smallest<Buffer - static_cast<int>(CountedAlready)>(varphi.get_saturation() + static_cast<int>(CountedAlready), count, HeapOrder{}));
        if(HeapOrder{}(disp, smallest)) {
          return PutInHeap::get(varphi, disp, get_type(point), get_type(other)) - PutInHeap::get(varphi, smallest, get_type(point), get_type(other));
        }
      }
    }
    return 0.;
  }

  template<int Buffer, bool CountedAlready, typename Point, typename Other, typename ToDiscard>
  static auto delta_discarding(const Saturated_model& varphi,
                               const ValueType& count,
                               const Point& point,
                               const Other& other,
                               const ToDiscard& to_discard) {
    // const auto sq_discard(normalized_square_distance(point, to_discard));
    // const auto discard_disp(PutInHeap::put(varphi, sq_discard, point, other));
    // ValueType count_copy(count);
    // count_copy.erase(std::remove(count_copy.begin(), count_copy.end(), discard_disp), count_copy.end());
    // std::make_heap(count_copy.begin(), count_copy.end(), HeapOrder{});
    // return delta<Buffer, CountedAlready>(varphi, count_copy, point, other);

    // TODO: Code below seems to work, but not super confident; test all functions with Catch!
    const auto sq(normalized_square_distance(point, other));
    if(sq >= varphi.get_square_lower_endpoint(get_type(point), get_type(other))) {
      if(count.size() < varphi.get_saturation() + static_cast<int>(CountedAlready)) {
        return varphi.apply(sq, get_type(point), get_type(other));
      } else if(count.size() == varphi.get_saturation() + static_cast<int>(CountedAlready)) {
        const auto sq_discard(normalized_square_distance(point, to_discard));
        const auto smallest(get_nth<0>(count, HeapOrder{}));
        if(sq_discard >= varphi.get_square_lower_endpoint(get_type(point), get_type(to_discard))) {
          const auto discard_disp(PutInHeap::put(varphi, sq_discard, point, to_discard));
          if(HeapOrder{}(discard_disp, smallest)) { // to_discard is in count
            return varphi.apply(sq, get_type(point), get_type(other));
          }
        }
        const auto disp(PutInHeap::put(varphi, sq, point, other));
        // to_discard is not in count
        if(HeapOrder{}(disp, smallest)) {
          return PutInHeap::get(varphi, disp, get_type(point), get_type(other)) - PutInHeap::get(varphi, smallest, get_type(point), get_type(other));
        }
      } else {
        const auto disp(PutInHeap::put(varphi, sq, point, other));
        const auto sq_discard(normalized_square_distance(point, to_discard));
        if(sq_discard >= varphi.get_square_lower_endpoint(get_type(point), get_type(to_discard))) {
          const auto discard_disp(PutInHeap::put(varphi, sq_discard, point, to_discard));
          const auto smallest(get_nth_smallest_excluding<Buffer - static_cast<int>(CountedAlready), 1>(varphi.get_saturation() + static_cast<int>(CountedAlready), count, discard_disp, HeapOrder{}));
          if(HeapOrder{}(disp, smallest)) {
            return PutInHeap::get(varphi, disp, get_type(point), get_type(other)) - PutInHeap::get(varphi, smallest, get_type(point), get_type(other));
          }
        } else {
          const auto smallest(get_nth_smallest<Buffer - static_cast<int>(CountedAlready)>(varphi.get_saturation() + static_cast<int>(CountedAlready), count, HeapOrder{}));
          if(HeapOrder{}(disp, smallest)) {
            return PutInHeap::get(varphi, disp, get_type(point), get_type(other)) - PutInHeap::get(varphi, smallest, get_type(point), get_type(other));
          }
        }
      }
    }
    return 0.;
  }
};

template<int Buffer, typename AbstractDispersion, int N = 1, typename Vector, typename CountVector, typename Point>
static void add_count_to_dispersion(const Saturated_model& varphi,
                                    Vector& dispersion,
                                    const CountVector& count_vector,
                                    const Point& point) {
  using size_t = typename Vector::size_type;
  using FloatType = typename Vector::value_type;
  for(size_t i(0); i < dispersion.size(); ++i) {
    const auto count_to_dispersion(static_cast<FloatType>(AbstractDispersion::template add_count_to_dispersion<Buffer>(varphi, count_vector, point, i)));
    dispersion[i] += static_cast<FloatType>(N) * count_to_dispersion;
  }
}

template<int Buffer, typename AbstractDispersion, int N = 1, typename Vector, typename CountVector, typename Point, typename ToDiscard>
static void add_count_to_dispersion_discarding(const Saturated_model& varphi,
                                               Vector& dispersion,
                                               const CountVector& count_vector,
                                               const Point& point,
                                               const ToDiscard& to_discard) {
  using size_t = typename Vector::size_type;
  using FloatType = typename Vector::value_type;
  for(size_t i(0); i < dispersion.size(); ++i) {
    FloatType count_to_dispersion{};
    if(i == static_cast<size_t>(get_type(to_discard))) {
      count_to_dispersion = static_cast<FloatType>(AbstractDispersion::template add_count_to_dispersion_discarding<Buffer>(varphi, count_vector, point, to_discard, i));
    } else {
      count_to_dispersion = static_cast<FloatType>(AbstractDispersion::template add_count_to_dispersion<Buffer>(varphi, count_vector, point, i));
    }
    dispersion[i] += static_cast<FloatType>(N) * count_to_dispersion;
  }
}

template<template<class> class T, typename... Args>
inline auto dispatch_model(const Saturated_model& model,
                           Args&&... args) {
  if(model.is_two_valued()) {
    return T<compute_dispersion_two_values_implementation>{}(model, std::forward<Args>(args)...);
  } else if(model.is_nonincreasing_after_lower_endpoint()) {
    return T<compute_dispersion_generic_implementation<put_square_distance>>{}(model, std::forward<Args>(args)...);
  } else {
    return T<compute_dispersion_generic_implementation<put_potential>>{}(model, std::forward<Args>(args)...);
  }
}

} // namespace detail
} // namespace ppjsdm

#endif // INCLUDE_COMPUTE_DISPERSION_IMPLEMENTATION
