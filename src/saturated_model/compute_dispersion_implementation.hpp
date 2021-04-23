#ifndef INCLUDE_COMPUTE_DISPERSION_IMPLEMENTATION
#define INCLUDE_COMPUTE_DISPERSION_IMPLEMENTATION

#include <Rcpp.h>

#include "saturated_model.hpp"
#include "../configuration/configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"
#include "../point/square_distance.hpp"
#include "../utility/flatten_strict_upper_triangular.hpp"

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

template<dispersionMethod Method>
struct compute_dispersion_implementation;

template<>
struct compute_dispersion_implementation<dispersionMethod::two_values> {
  using ValueType = unsigned long long int;

  template<typename Vector, typename Point>
  static auto add_count_to_dispersion(const Saturated_model& varphi,
                                      const Vector& count_vector,
                                      const Point&,
                                      decltype(count_vector.size()) i) {
    return std::min(varphi.get_saturation(), count_vector[i]);
  }

  template<typename Point, typename ToDiscard, typename Vector>
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

  template<typename Point, typename Other>
  static void update_count_nonsaturated(const Saturated_model& varphi, ValueType& count, const Point& point, const Other& other) {
    if(apply_potential(varphi, point, other) > 0.) {
      ++count;
    }
  }

  template<bool CountedAlready, typename Point, typename Other>
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

  template<bool CountedAlready, typename Point, typename Other, typename ToDiscard>
  static auto delta_discarding(const Saturated_model& varphi,
                               const ValueType& count,
                               const Point& point,
                               const Other& other,
                               const ToDiscard& to_discard) {
    if((count < varphi.get_saturation() + static_cast<int>(CountedAlready) + static_cast<int>(apply_potential(varphi, point, to_discard) > 0.)) && (apply_potential(varphi, point, other) > 0.)) {
      return 1.;
    } else {
      return 0.;
    }
  }
};

template<>
struct compute_dispersion_implementation<dispersionMethod::nonincreasing_after_lower_endpoint> {
  using ValueType = std::vector<double>;

  template<typename Point, typename Vector>
  static auto add_count_to_dispersion(const Saturated_model& varphi,
                                      const Vector& count_vector,
                                      const Point& point,
                                      decltype(count_vector.size()) i) {
    if(count_vector[i].size() <= varphi.get_saturation()) {
      return std::accumulate(count_vector[i].begin(),
                             count_vector[i].end(),
                             0., [&varphi, i, &point](double count, auto val) {
                               return count + varphi.apply(val, i, get_type(point));
                             });
    } else if(count_vector[i].size() == varphi.get_saturation() + 1) {
      return std::accumulate(count_vector[i].begin() + 1,
                             count_vector[i].end(),
                             0., [&varphi, i, &point](double count, auto val) {
                               return count + varphi.apply(val, i, get_type(point));
                             });
    } else {
      const auto to_skip(std::max(count_vector[i][1], count_vector[i][2]));
      return std::accumulate(count_vector[i].begin() + 1,
                             count_vector[i].end(),
                             0., [&varphi, i, &point, to_skip](double count, auto val) {
                               if(val != to_skip) {
                                 return count + varphi.apply(val, i, get_type(point));
                               } else {
                                 return count;
                               }
                             });
    }
  }

  template<typename Point, typename ToDiscard, typename Vector>
  static auto add_count_to_dispersion_discarding(const Saturated_model& varphi,
                                                 const Vector& count_vector,
                                                 const Point& point,
                                                 const ToDiscard& to_discard,
                                                 decltype(count_vector.size()) i) {
    // TODO: Temp solution to make sure everything works, it's extremely inefficient
    const auto discard_sq(normalized_square_distance(point, to_discard));
    Vector count_vector_copy(count_vector);
    count_vector_copy[i].erase(std::remove(count_vector_copy[i].begin(), count_vector_copy[i].end(), discard_sq), count_vector_copy[i].end());
    std::make_heap(count_vector_copy[i].begin(), count_vector_copy[i].end());

    return add_count_to_dispersion(varphi, count_vector_copy, point, i);
  }

  template<int Buffer = 1, typename Point, typename Other>
  static void update_count(const Saturated_model& varphi, ValueType& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(sq >= varphi.get_square_lower_endpoint(get_type(point), get_type(other))) {
      if(count.size() < varphi.get_saturation() + Buffer) {
        count.emplace_back(sq);
        std::push_heap(count.begin(), count.end());
      } else if(sq < count[0]) {
        count.emplace_back(sq);
        std::pop_heap(count.begin(), count.end());
        count.pop_back();
      }
    }
  }

  template<typename Point, typename Other>
  static void update_count_nonsaturated(const Saturated_model& varphi, ValueType& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(sq >= varphi.get_square_lower_endpoint(get_type(point), get_type(other))) {
      count.emplace_back(sq);
      std::push_heap(count.begin(), count.end());
    }
  }

  template<bool CountedAlready, typename Point, typename Other>
  static auto delta(const Saturated_model& varphi,
                    const ValueType& count,
                    const Point& point,
                    const Other& other) {
    // const auto sq(normalized_square_distance(point, other));
    // if(sq >= varphi.get_square_lower_endpoint(get_type(point), get_type(other))) {
    //   // TODO: if constexpr
    //   if(CountedAlready) {
    //     if(count.size() <= varphi.get_saturation()) {
    //       dispersion += varphi.apply(sq, get_type(point), get_type(other));
    //     } else if(sq < count[0]) {
    //       dispersion += varphi.apply(sq, get_type(point), get_type(other)) - varphi.apply(count[0], get_type(point), get_type(other));
    //     }
    //   } else {
    //     if(count.size() < varphi.get_saturation()) {
    //       dispersion += varphi.apply(sq, get_type(point), get_type(other));
    //     } else {
    //       const auto m(count.size() == varphi.get_saturation() ? count[0] : std::max(count[1], count[2]));
    //       if(sq < m) {
    //         dispersion += varphi.apply(sq, get_type(point), get_type(other)) - varphi.apply(m, get_type(point), get_type(other));
    //       }
    //     }
    //   }
    // }
    const auto sq(normalized_square_distance(point, other));
    if(sq >= varphi.get_square_lower_endpoint(get_type(point), get_type(other))) {
      // TODO: if constexpr
      if(CountedAlready) {
        if(count.size() <= varphi.get_saturation()) {
          return varphi.apply(sq, get_type(point), get_type(other));
        } else {
          const auto largest(count.size() == varphi.get_saturation() + 1 ? count[0] : std::max(count[std::min<int>(1, count.size() - 1)], count[std::min<int>(2, count.size() - 1)]));
          if(sq < largest) {
            return varphi.apply(sq, get_type(point), get_type(other)) - varphi.apply(largest, get_type(point), get_type(other));
          }
        }
      } else {
        if(count.size() < varphi.get_saturation()) {
          return varphi.apply(sq, get_type(point), get_type(other));
        } else {
          // TODO: Horrible code below
          const auto largest(count.size() == varphi.get_saturation() ?
                               count[0] :
                               count.size() == varphi.get_saturation() + 1 ?
                               std::max(count[std::min<int>(1, count.size() - 1)], count[std::min<int>(2, count.size() - 1)]) :
                               std::max(std::min(count[std::min<int>(1, count.size() - 1)], count[std::min<int>(2, count.size() - 1)]),
                                        std::max(std::max(count[std::min<int>(3, count.size() - 1)], count[std::min<int>(4, count.size() - 1)]),
                                                 std::max(count[std::min<int>(5, count.size() - 1)], count[std::min<int>(6, count.size() - 1)]))));
          if(sq < largest) {
            return varphi.apply(sq, get_type(point), get_type(other)) - varphi.apply(largest, get_type(point), get_type(other));
          }
        }
      }
    }
    return 0.;
  }

  template<bool CountedAlready, typename Point, typename Other, typename ToDiscard>
  static auto delta_discarding(const Saturated_model& varphi,
                               const ValueType& count,
                               const Point& point,
                               const Other& other,
                               const ToDiscard& to_discard) {
    const auto discard_sq(normalized_square_distance(point, to_discard));
    // TODO: Temp solution to make sure everything works, it's extremely inefficient
    ValueType count_copy(count);
    count_copy.erase(std::remove(count_copy.begin(), count_copy.end(), discard_sq), count_copy.end());
    std::make_heap(count_copy.begin(), count_copy.end());

    return delta<CountedAlready>(varphi, count_copy, point, other);
  }
};

template<>
struct compute_dispersion_implementation<dispersionMethod::generic> {
  using ValueType = std::vector<double>;

  template<typename Point, typename Vector>
  static auto add_count_to_dispersion(const Saturated_model& varphi,
                                      const Vector& count_vector,
                                      const Point&,
                                      decltype(count_vector.size()) i) {
    if(count_vector[i].size() <= varphi.get_saturation()) {
      return std::accumulate(count_vector[i].begin(),
                             count_vector[i].end(),
                             0.);
    } else if(count_vector[i].size() == varphi.get_saturation() + 1) {
      return std::accumulate(count_vector[i].begin() + 1,
                             count_vector[i].end(),
                             0.);
    } else {
      const auto to_skip(std::min(count_vector[i][1], count_vector[i][2]));
      return std::accumulate(count_vector[i].begin() + 1,
                             count_vector[i].end(),
                             0., [to_skip](double count, auto val) {
                               if(val != to_skip) {
                                 return count + val;
                               } else {
                                 return count;
                               }
                             });
    }
  }

  template<typename Point, typename ToDiscard, typename Vector>
  static auto add_count_to_dispersion_discarding(const Saturated_model& varphi,
                                                 const Vector& count_vector,
                                                 const Point& point,
                                                 const ToDiscard& to_discard,
                                                 decltype(count_vector.size()) i) {
    // TODO: Temp solution to make sure everything works, it's extremely inefficient
    const auto discard_disp(apply_potential(varphi, point, to_discard));
    Vector count_vector_copy(count_vector);
    count_vector_copy[i].erase(std::remove(count_vector_copy[i].begin(), count_vector_copy[i].end(), discard_disp), count_vector_copy[i].end());
    std::make_heap(count_vector_copy[i].begin(), count_vector_copy[i].end(), std::greater<double>{});

    return add_count_to_dispersion(varphi, count_vector_copy, point, i);
  }

  template<int Buffer = 1, typename Point, typename Other>
  static void update_count(const Saturated_model& varphi, ValueType& count, const Point& point, const Other& other) {
    const auto disp(apply_potential(varphi, other, point));
    if(count.size() < varphi.get_saturation() + Buffer) {
      count.emplace_back(disp);
      std::push_heap(count.begin(), count.end(), std::greater<double>{});
    } else if(disp > count[0]) {
      count.emplace_back(disp);
      std::pop_heap(count.begin(), count.end(), std::greater<double>{});
      count.pop_back();
    }
  }

  template<typename Point, typename Other>
  static void update_count_nonsaturated(const Saturated_model& varphi, ValueType& count, const Point& point, const Other& other) {
    const auto disp(apply_potential(varphi, other, point));
    count.emplace_back(disp);
    std::push_heap(count.begin(), count.end(), std::greater<double>{});
  }

  template<bool CountedAlready, typename Point, typename Other>
  static auto delta(const Saturated_model& varphi,
                    const ValueType& count,
                    const Point& point,
                    const Other& other) {
    // const auto disp(apply_potential(varphi, other, point));
    // // TODO: if constexpr
    // if(CountedAlready) {
    //   if(count.size() <= varphi.get_saturation()) {
    //     dispersion += disp;
    //   } else if(disp > count[0]) {
    //     dispersion += disp - count[0];
    //   }
    // } else {
    //   if(count.size() < varphi.get_saturation()) {
    //     dispersion += disp;
    //   } else {
    //     const auto m(count.size() == varphi.get_saturation() ? count[0] : std::min(count[1], count[2]));
    //     if(disp > m) {
    //       dispersion += disp - m;
    //     }
    //   }
    // }
    const auto disp(apply_potential(varphi, other, point));
    // TODO: if constexpr
    if(CountedAlready) {
      if(count.size() <= varphi.get_saturation()) {
        return disp;
      } else {
        const auto smallest(count.size() == varphi.get_saturation() + 1 ? count[0] : std::min(count[std::min<int>(1, count.size() - 1)], count[std::min<int>(2, count.size() - 1)]));
        if(disp > smallest) {
          return disp - smallest;
        }
      }
    } else {
      if(count.size() < varphi.get_saturation()) {
        return disp;
      } else {
        // TODO: Horrible code below
        const auto smallest(count.size() == varphi.get_saturation() ?
                             count[0] :
                             count.size() == varphi.get_saturation() + 1 ?
                             std::min(count[std::min<int>(1, count.size() - 1)], count[std::min<int>(2, count.size() - 1)]) :
                             std::min(std::max(count[std::min<int>(1, count.size() - 1)], count[std::min<int>(2, count.size() - 1)]),
                                      std::min(std::min(count[std::min<int>(3, count.size() - 1)], count[std::min<int>(4, count.size() - 1)]),
                                               std::min(count[std::min<int>(5, count.size() - 1)], count[std::min<int>(6, count.size() - 1)]))));
        if(disp > smallest) {
          return disp - smallest;
        }
      }
    }
    return 0.;
  }

  template<bool CountedAlready, typename Point, typename Other, typename ToDiscard>
  static auto delta_discarding(const Saturated_model& varphi,
                               const ValueType& count,
                               const Point& point,
                               const Other& other,
                               const ToDiscard& to_discard) {
    const auto discard_disp(apply_potential(varphi, point, to_discard));
    // TODO: Temp solution to make sure everything works, it's extremely inefficient
    ValueType count_copy(count);
    count_copy.erase(std::remove(count_copy.begin(), count_copy.end(), discard_disp), count_copy.end());
    std::make_heap(count_copy.begin(), count_copy.end(), std::greater<double>{});

    return delta<CountedAlready>(varphi, count_copy, point, other);
  //   const auto disp(apply_potential(varphi, other, point));
  //   const auto disp_to_discard(apply_potential(varphi, point, to_discard));
  //   const bool is_discard_in_count(std::find(count.begin(), count.end(), disp_to_discard) != count.end());
  //   // TODO: if constexpr
  //   if(CountedAlready) {
  //     if(is_discard_in_count) {
  //       if(count.size() <= 1 + varphi.get_saturation()) {
  //         dispersion += disp;
  //       } else {
  //         const auto smallest(std::min(disp_to_discard, count[0]));
  //         if(disp > smallest) {
  //           dispersion += disp - smallest;
  //         }
  //       }
  //     } else {
  //       if(count.size() <= varphi.get_saturation()) {
  //         dispersion += disp;
  //       } else if(disp > count[0]) {
  //         dispersion += disp - count[0];
  //       }
  //     }
  //   } else {
  //     if(is_discard_in_count) {
  //       if(count.size() < 1 + varphi.get_saturation()) {
  //         dispersion += disp;
  //       } else {
  //         auto smallest(count[0]);
  //         auto second_smallest(std::min(count[1], count[2]));
  //         if(disp_to_discard == smallest) {
  //           smallest = second_smallest;
  //           second_smallest = std::min(std::min(count[3], count[4]), std::min(count[5], count[6]));
  //         } else if(disp_to_discard == second_smallest) {
  //           second_smallest = std::min(std::min(count[3], count[4]), std::min(count[5], count[6]));
  //         }
  //         const auto m(count.size() == varphi.get_saturation() ? smallest : second_smallest);
  //         if(disp > m) {
  //           dispersion += disp - m;
  //         }
  //       }
  //     } else {
  //       if(count.size() < varphi.get_saturation()) {
  //         dispersion += disp;
  //       } else {
  //         const auto m(count.size() == varphi.get_saturation() ? count[0] : std::min(count[1], count[2]));
  //         if(disp > m) {
  //           dispersion += disp - m;
  //         }
  //       }
  //     }
  //   }
  }
};

template<typename AbstractDispersion, int N = 1, typename Vector, typename CountVector, typename Point>
static void add_count_to_dispersion(const Saturated_model& varphi,
                                    Vector& dispersion,
                                    const CountVector& count_vector,
                                    const Point& point) {
  using size_t = typename Vector::size_type;
  using FloatType = typename Vector::value_type;
  for(size_t i(0); i < dispersion.size(); ++i) {
    const auto count_to_dispersion(static_cast<FloatType>(AbstractDispersion::add_count_to_dispersion(varphi, count_vector, point, i)));
    dispersion[i] += static_cast<FloatType>(N) * count_to_dispersion;
  }
}

template<typename AbstractDispersion, int N = 1, typename Vector, typename CountVector, typename Point, typename ToDiscard>
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
      count_to_dispersion = static_cast<FloatType>(AbstractDispersion::add_count_to_dispersion_discarding(varphi, count_vector, point, to_discard, i));
    } else {
      count_to_dispersion = static_cast<FloatType>(AbstractDispersion::add_count_to_dispersion(varphi, count_vector, point, i));
    }
    dispersion[i] += static_cast<FloatType>(N) * count_to_dispersion;
  }
}

} // namespace detail
} // namespace ppjsdm

#endif // INCLUDE_COMPUTE_DISPERSION_IMPLEMENTATION
