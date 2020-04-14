#ifndef INCLUDE_PPJSDM_SATURATED_MODEL
#define INCLUDE_PPJSDM_SATURATED_MODEL

#include <Rcpp.h>

#include "potentials/medium_range_potentials.hpp"
#include "potentials/short_range_potentials.hpp"
#include "../configuration/configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"
#include "../point/square_distance.hpp"
#include "../utility/for_each_container.hpp"

#include <algorithm> // std::accumulate, std::pop_heap, std::push_heap
#include <functional> // std::greater
#include <utility> // std::forward, std::move
#include <vector> // std::vector

namespace ppjsdm {
/*namespace detail {

template<typename Point, typename Dispersion, typename Select = void>
class Compute_dispersion_implementation;

template<typename Point, typename Dispersion>
class Compute_dispersion_implementation<Point, Dispersion, std::enable_if_t<has_nonzero_value_v<Dispersion>>> {
protected:
  Compute_dispersion_implementation(const Point& point,
                                    R_xlen_t number_types,
                                    const Dispersion& dispersion):
  point_(point),
  dispersion_(dispersion),
  count_vector_(number_types) {}

  using CountType = std::vector<unsigned long long int>;
  const Point& point_;
  const Dispersion& dispersion_;
  CountType count_vector_;

  template<typename Configuration>
  void update_count(const Configuration& configuration) {
    for(const auto& current_point: configuration) {
      if(count_vector_[get_type(current_point)] < dispersion_.saturation_ && dispersion_.apply(current_point, point_) > 0) {
        ++count_vector_[get_type(current_point)];
      }
    }
  }

  auto compute_dispersion(const CountType& count) const {
    const auto number_types(count.size());
    std::vector<double> dispersion(number_types);
    using size_t = typename decltype(dispersion)::size_type;
    constexpr double nonzero_value(Dispersion::nonzero_value);
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      dispersion[i] = nonzero_value * static_cast<double>(count[i]);
    }
    return dispersion;
  }
};

template<typename Point, typename Dispersion>
class Compute_dispersion_implementation<Point, Dispersion,
                                        std::enable_if_t<!has_nonzero_value_v<Dispersion>
&& Dispersion::is_nonincreasing>> {
protected:
  Compute_dispersion_implementation(const Point& point,
                                    R_xlen_t number_types,
                                    const Dispersion& dispersion):
  point_(point),
  dispersion_(dispersion),
  count_vector_(number_types) {}

  using CountType = std::vector<std::vector<double>>;
  const Point& point_;
  const Dispersion& dispersion_;
  CountType count_vector_;

  template<typename Configuration>
  void update_count(const Configuration& configuration) {
    for(const auto& current_point: configuration) {
      const auto sq(normalized_square_distance(current_point, point_));
      auto& current(count_vector_[get_type(current_point)]);
      if(current.size() < dispersion_.saturation_) {
        current.emplace_back(sq);
        std::push_heap(current.begin(), current.end());
      } else if(sq < current[0]) {
        current.emplace_back(sq);
        std::pop_heap(current.begin(), current.end());
        current.pop_back();
      }
    }
  }

  auto compute_dispersion(const CountType& count) const {
    const auto number_types(count.size());
    std::vector<double> dispersion(number_types);
    using size_t = typename decltype(dispersion)::size_type;
    const auto point_type(get_type(point_));
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      double d(0);
      for(const auto& c: count[i]) {
        d += dispersion_.apply(c, i, point_type);
      }
      dispersion[i] = d;
    }
    return dispersion;
  }
};

template<typename Point, typename Dispersion>
class Compute_dispersion_implementation<Point, Dispersion,
                                        std::enable_if_t<!has_nonzero_value_v<Dispersion>
&& !Dispersion::is_nonincreasing
&& is_nonincreasing_after_lower_endpoint_v<Dispersion>>> {
protected:
  Compute_dispersion_implementation(const Point& point,
                                    R_xlen_t number_types,
                                    const Dispersion& dispersion):
  point_(point),
  dispersion_(dispersion),
  count_vector_(number_types) {}

  using CountType = std::vector<std::vector<double>>;
  const Point& point_;
  const Dispersion& dispersion_;
  CountType count_vector_;

  template<typename Configuration>
  void update_count(const Configuration& configuration) {
    for(const auto& current_point: configuration) {
      const auto sq(normalized_square_distance(current_point, point_));
      // TODO: This doesn't really change from the case above; factorise?
      if(sq >= Dispersion::get_square_lower_endpoint(get_type(current_point), get_type(point_))) {
        auto& current(count_vector_[get_type(current_point)]);
        if(current.size() < dispersion_.saturation_) {
          current.emplace_back(sq);
          std::push_heap(current.begin(), current.end());
        } else if(sq < current[0]) {
          current.emplace_back(sq);
          std::pop_heap(current.begin(), current.end());
          current.pop_back();
        }
      }
    }
  }

  auto compute_dispersion(const CountType& count) const {
    const auto number_types(count.size());
    std::vector<double> dispersion(number_types);
    using size_t = typename decltype(dispersion)::size_type;
    const auto point_type(get_type(point_));
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      double d(0);
      for(const auto& c: count[i]) {
        d += dispersion_.apply(c, i, point_type);
      }
      dispersion[i] = d;
    }
    return dispersion;
  }
};

template<typename Point, typename Dispersion>
class Compute_dispersion_implementation<Point, Dispersion,
                                        std::enable_if_t<!has_nonzero_value_v<Dispersion>
&& !Dispersion::is_nonincreasing
&& !is_nonincreasing_after_lower_endpoint_v<Dispersion>>> {
protected:
  Compute_dispersion_implementation(const Point& point,
                                    R_xlen_t number_types,
                                    const Dispersion& dispersion):
  point_(point),
  dispersion_(dispersion),
  count_vector_(number_types) {}

  using CountType = std::vector<std::vector<double>>;
  const Point& point_;
  const Dispersion& dispersion_;
  CountType count_vector_;

  template<typename Configuration>
  void update_count(const Configuration& configuration) {
    for(const auto& current_point: configuration) {
      const auto disp(dispersion_.apply(current_point, point_));
      auto& current(count_vector_[get_type(current_point)]);
      if(current.size() < dispersion_.saturation_) {
        current.emplace_back(disp);
        std::push_heap(current.begin(), current.end(), std::greater<double>{});
      } else if(disp > current[0]) {
        current.emplace_back(disp);
        std::pop_heap(current.begin(), current.end(), std::greater<double>{});
        current.pop_back();
      }
    }
  }

  auto compute_dispersion(const CountType& count) const {
    const auto number_types(count.size());
    std::vector<double> dispersion(number_types);
    using size_t = typename decltype(dispersion)::size_type;
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      dispersion[i] = std::accumulate(count[i].begin(), count[i].end(), 0.);
    }
    return dispersion;
  }
};

} // namespace detail

template<typename Point, typename Dispersion>
class Compute_dispersion : public detail::Compute_dispersion_implementation<Point, Dispersion> {
private:
  using ImplementationBase = detail::Compute_dispersion_implementation<Point, Dispersion>;
public:
  template<typename... Args>
  Compute_dispersion(Args&&... args):
    ImplementationBase(std::forward<Args>(args)...) {}

  template<typename Configuration>
  void add_configuration(const Configuration& configuration) {
    ImplementationBase::update_count(configuration);
  }

  auto compute() const {
    return ImplementationBase::compute_dispersion(ImplementationBase::count_vector_);
  }

  auto compute_from_state(const typename ImplementationBase::CountType& count) const {
    return ImplementationBase::compute_dispersion(count);
  }

  const auto& get_state() const {
    return ImplementationBase::count_vector_;
  }
};

// Note: Use inheritance to benefit from EBO.
template<typename Varphi>
class Saturated_model: public Varphi {
private:
  template<typename, typename, typename>
  friend class detail::Compute_dispersion_implementation;
public:
  template<typename... Args>
  Saturated_model(unsigned long long int saturation, Args&&... args):
    Varphi(std::forward<Args>(args)...), saturation_(saturation) {}

  template<typename Point>
  auto get_compute_dispersion_object(const Point& point,
                                     R_xlen_t number_types) const {
    return Compute_dispersion<Point, Saturated_model<Varphi>>(point, number_types, *this);
  }

  template<typename Point, typename Configuration>
  auto compute(const Point& point,
               R_xlen_t number_types,
               const Configuration& configuration) const {
    auto dispersion(get_compute_dispersion_object(point, number_types));
    dispersion.add_configuration(configuration);
    return dispersion.compute();
  }

  template<typename Window>
  double get_maximum(const Window& window) const {
    return static_cast<double>(saturation_) * Varphi::get_maximum(window);
  }
private:
  unsigned long long int saturation_;
};

namespace detail {

template<typename Point, typename Dispersion, typename Select = void>
class Fixed_compute_dispersion_implementation;

template<typename Point, typename Dispersion>
class Fixed_compute_dispersion_implementation<Point, Dispersion, std::enable_if_t<has_nonzero_value_v<Dispersion>>> {
protected:
  Fixed_compute_dispersion_implementation(const Point& point,
                                    R_xlen_t number_types,
                                    const Dispersion& dispersion):
  point_(point),
  dispersion_(dispersion),
  count_vector_(number_types) {}

  using CountType = std::vector<unsigned long long int>;
  const Point& point_;
  const Dispersion& dispersion_;
  CountType count_vector_;

  template<typename Configuration>
  void update_count(const Configuration& configuration) {
    for(const auto& current_point: configuration) {
      if(count_vector_[get_type(current_point)] < dispersion_.saturation_ && dispersion_.apply(current_point, point_) > 0) {
        ++count_vector_[get_type(current_point)];
      }
    }
  }

  auto compute_dispersion(const CountType& count) const {
    const auto number_types(count.size());
    std::vector<double> dispersion(number_types);
    using size_t = typename decltype(dispersion)::size_type;
    constexpr double nonzero_value(Dispersion::nonzero_value);
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      dispersion[i] = nonzero_value * static_cast<double>(count[i]);
    }
    return dispersion;
  }
};

template<typename Point, typename Dispersion>
class Fixed_compute_dispersion_implementation<Point, Dispersion,
                                        std::enable_if_t<!has_nonzero_value_v<Dispersion>
&& Dispersion::is_nonincreasing>> {
protected:
  Fixed_compute_dispersion_implementation(const Point& point,
                                    R_xlen_t number_types,
                                    const Dispersion& dispersion):
  point_(point),
  dispersion_(dispersion),
  count_vector_(number_types) {}

  using CountType = std::vector<std::vector<double>>;
  const Point& point_;
  const Dispersion& dispersion_;
  CountType count_vector_;

  template<typename Configuration>
  void update_count(const Configuration& configuration) {
    for(const auto& current_point: configuration) {
      const auto sq(normalized_square_distance(current_point, point_));
      auto& current(count_vector_[get_type(current_point)]);
      if(current.size() < dispersion_.saturation_) {
        current.emplace_back(sq);
        std::push_heap(current.begin(), current.end());
      } else if(sq < current[0]) {
        current.emplace_back(sq);
        std::pop_heap(current.begin(), current.end());
        current.pop_back();
      }
    }
  }

  auto compute_dispersion(const CountType& count) const {
    const auto number_types(count.size());
    std::vector<double> dispersion(number_types);
    using size_t = typename decltype(dispersion)::size_type;
    const auto point_type(get_type(point_));
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      double d(0);
      for(const auto& c: count[i]) {
        d += dispersion_.apply(c, i, point_type);
      }
      dispersion[i] = d;
    }
    return dispersion;
  }
};

template<typename Point, typename Dispersion>
class Fixed_compute_dispersion_implementation<Point, Dispersion,
                                        std::enable_if_t<!has_nonzero_value_v<Dispersion>
&& !Dispersion::is_nonincreasing
&& is_nonincreasing_after_lower_endpoint_v<Dispersion>>> {
protected:
  Fixed_compute_dispersion_implementation(const Point& point,
                                    R_xlen_t number_types,
                                    const Dispersion& dispersion):
  point_(point),
  dispersion_(dispersion),
  count_vector_(number_types) {}

  using CountType = std::vector<std::vector<double>>;
  const Point& point_;
  const Dispersion& dispersion_;
  CountType count_vector_;

  template<typename Configuration>
  void update_count(const Configuration& configuration) {
    for(const auto& current_point: configuration) {
      const auto sq(normalized_square_distance(current_point, point_));
      // TODO: This doesn't really change from the case above; factorise?
      if(sq >= Dispersion::get_square_lower_endpoint(get_type(current_point), get_type(point_))) {
        auto& current(count_vector_[get_type(current_point)]);
        if(current.size() < dispersion_.saturation_) {
          current.emplace_back(sq);
          std::push_heap(current.begin(), current.end());
        } else if(sq < current[0]) {
          current.emplace_back(sq);
          std::pop_heap(current.begin(), current.end());
          current.pop_back();
        }
      }
    }
  }

  auto compute_dispersion(const CountType& count) const {
    const auto number_types(count.size());
    std::vector<double> dispersion(number_types);
    using size_t = typename decltype(dispersion)::size_type;
    const auto point_type(get_type(point_));
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      double d(0);
      for(const auto& c: count[i]) {
        d += dispersion_.apply(c, i, point_type);
      }
      dispersion[i] = d;
    }
    return dispersion;
  }
};

template<typename Point, typename Dispersion>
class Fixed_compute_dispersion_implementation<Point, Dispersion,
                                        std::enable_if_t<!has_nonzero_value_v<Dispersion>
&& !Dispersion::is_nonincreasing
&& !is_nonincreasing_after_lower_endpoint_v<Dispersion>>> {
protected:
  Fixed_compute_dispersion_implementation(const Point& point,
                                    R_xlen_t number_types,
                                    const Dispersion& dispersion):
  point_(point),
  dispersion_(dispersion),
  count_vector_(number_types) {}

  using CountType = std::vector<std::vector<double>>;
  const Point& point_;
  const Dispersion& dispersion_;
  CountType count_vector_;

  template<typename Configuration>
  void update_count(const Configuration& configuration) {
    for(const auto& current_point: configuration) {
      const auto disp(dispersion_.apply(current_point, point_));
      auto& current(count_vector_[get_type(current_point)]);
      if(current.size() < dispersion_.saturation_) {
        current.emplace_back(disp);
        std::push_heap(current.begin(), current.end(), std::greater<double>{});
      } else if(disp > current[0]) {
        current.emplace_back(disp);
        std::pop_heap(current.begin(), current.end(), std::greater<double>{});
        current.pop_back();
      }
    }
  }

  auto compute_dispersion(const CountType& count) const {
    const auto number_types(count.size());
    std::vector<double> dispersion(number_types);
    using size_t = typename decltype(dispersion)::size_type;
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      dispersion[i] = std::accumulate(count[i].begin(), count[i].end(), 0.);
    }
    return dispersion;
  }
};

} // namespace detail

template<typename Point, typename Dispersion>
class Fixed_compute_dispersion : public detail::Fixed_compute_dispersion_implementation<Point, Dispersion> {
private:
  using ImplementationBase = detail::Fixed_compute_dispersion_implementation<Point, Dispersion>;
public:
  template<typename... Args>
  Fixed_compute_dispersion(Args&&... args):
    ImplementationBase(std::forward<Args>(args)...) {}

  template<typename Configuration>
  void add_configuration(const Configuration& configuration) {
    ImplementationBase::update_count(configuration);
  }

  auto compute_from_state(const typename ImplementationBase::CountType& count) const {
    return ImplementationBase::compute_dispersion(count);
  }

  auto compute() const {
    return compute_from_state(ImplementationBase::count_vector_);
  }

  const auto& get_state() const {
    return ImplementationBase::count_vector_;
  }
};

// Note: Use inheritance to benefit from EBO.
template<typename Varphi>
class Fixed_saturated_model: public Varphi {
private:
  template<typename, typename, typename>
  friend class detail::Fixed_compute_dispersion_implementation;
public:
  template<typename... Args>
  Fixed_saturated_model(unsigned long long int saturation, Args&&... args):
    Varphi(std::forward<Args>(args)...),
    saturation_(saturation) {}

  template<typename Point>
  auto get_compute_dispersion_object(const Point& point,
                                     R_xlen_t number_types) const {
    return Fixed_compute_dispersion<Point, Fixed_saturated_model<Varphi>>(point, number_types, *this);
  }

  template<typename Point, typename Configuration>
  auto compute(const Point& point,
               R_xlen_t number_types,
               const Configuration& configuration) const {
    auto dispersion(get_compute_dispersion_object(point, number_types));
    dispersion.add_configuration(configuration);
    return dispersion.compute();
  }

  template<typename Window>
  double get_maximum(const Window& window) const {
    return Varphi::get_maximum(window) * static_cast<double>(saturation_);
  }
private:
  unsigned long long int saturation_;
};

// TODO: When everything is fixed, note that this add_configuration/compute_from_state system isn't useful.
template<typename Point>
class Compute_Dixon_dispersion {
public:
  Compute_Dixon_dispersion(const Point& point, R_xlen_t number_types, unsigned long long int saturation):
  point_(point), number_types_(number_types), points_{}, saturation_(saturation) {}

  template<typename Configuration>
  void add_configuration(const Configuration& configuration) {
    points_.reserve(size(configuration));
    for(const auto& point: configuration) {
      add_point(points_, point);
    }
  }

  auto compute() const {
    return compute(points_);
  }

  auto compute_from_state(const std::vector<Marked_point>& points) const {
    return compute(points);
  }

  const auto& get_state() const {
    return points_;
  }
private:
  Point point_;
  R_xlen_t number_types_;
  std::vector<Marked_point> points_;
  unsigned long long int saturation_;

  auto compute(std::vector<Marked_point> points) const {
    const auto segregation(compute_segregation(points, get_type(point_)));
    points.push_back(point_);
    auto segregation_plus(compute_segregation(points, get_type(point_)));
    std::transform(segregation_plus.begin(), segregation_plus.end(), segregation.begin(),
                   segregation_plus.begin(), std::minus<double>());
    return segregation_plus;
  }

  auto compute_segregation(const std::vector<Marked_point>& points, int type) const {
    std::vector<unsigned long long int> weights(number_types_);
    std::vector<unsigned long long int> closest_types(number_types_);
    for(const auto& current_point: points) {
      const auto current_type(get_type(current_point));
      ++weights[current_type];
      if(current_type == type) {
        double closest_square_distance(std::numeric_limits<double>::infinity());
        auto closest_type(current_type);
        for(const auto& other_point: points) {
          if(other_point != current_point) {
            const auto sq(square_distance(current_point, other_point));
            if(closest_square_distance > sq) {
              closest_square_distance = sq;
              closest_type = get_type(other_point);
            }
          }
        }
        ++closest_types[closest_type];
      }
    }
    std::vector<double> segregation(number_types_);
    for(R_xlen_t i(0); i < number_types_; ++i) {
      if(i == type) {
        if(closest_types[i] > 0U && weights[i] > closest_types[i] && static_cast<unsigned long long int>(size(points)) > weights[i]) {
          const auto value(std::log(closest_types[i] / (weights[i] - closest_types[i])) / ((weights[i] - 1) / (size(points) - weights[i])));
          const auto truncated(std::max(std::min(value, static_cast<double>(saturation_) / 2.0), -static_cast<double>(saturation_) / 2.0));
          segregation[i] = truncated + static_cast<double>(saturation_) / 2.0;
        } else {
          segregation[i] = static_cast<double>(saturation_) / 2.0;
        }
      } else {
        if(closest_types[i] > 0U && weights[type] > closest_types[i] && static_cast<unsigned long long int>(size(points)) > weights[i] + 1U) {
          const auto value(std::log(closest_types[i] / (weights[type] - closest_types[i])) / (weights[type] / (size(points) - weights[i] - 1)));
          const auto truncated(std::max(std::min(value, static_cast<double>(saturation_) / 2.0), -static_cast<double>(saturation_) / 2.0));
          segregation[i] = truncated + static_cast<double>(saturation_) / 2.0;
        } else {
          segregation[i] = static_cast<double>(saturation_) / 2.0;
        }
      }
    }
    return segregation;
  }
};

class Dixon_model {
public:
  Dixon_model(unsigned long long int saturation):
  saturation_(saturation) {}

  template<typename Point>
  auto get_compute_dispersion_object(const Point& point,
                                     R_xlen_t number_types) const {
    return Compute_Dixon_dispersion<Point>(point, number_types, saturation_);
  }

  template<typename Point, typename Configuration>
  auto compute(const Point& point,
               R_xlen_t number_types,
               const Configuration& configuration) const {
    auto dispersion(get_compute_dispersion_object(point, number_types));
    dispersion.add_configuration(configuration);
    return dispersion.compute();
  }

  template<typename Window>
  double get_maximum(const Window&) const {
    return static_cast<double>(saturation_);
  }
private:
  unsigned long long int saturation_;
};*/

class Saturated_model {
public:
  virtual bool is_nonincreasing() const = 0;
  virtual bool is_nonincreasing_after_lower_endpoint() const = 0;
  virtual bool is_two_valued() const = 0;
  virtual double apply(double normalized_square_distance, int i, int j) const = 0;
  virtual double get_maximum() const = 0;
  virtual unsigned long long int get_saturation() const = 0;
  virtual double get_square_lower_endpoint(int i, int j) const = 0;
};

template<typename Point, typename Other>
inline auto apply_potential(const Saturated_model& potential, const Point& point, const Other& other) {
  return potential.apply(normalized_square_distance(point, other), get_type(point), get_type(other));
}

namespace detail {

enum class dispersionMethod {two_values, nonincreasing, nonincreasing_after_lower_endpoint, generic};

template<bool Approximate, dispersionMethod Method>
class compute_dispersion;

template<bool Approximate>
class compute_dispersion<Approximate, dispersionMethod::two_values> {
public:
  compute_dispersion() {}

  template<typename Point, typename... Configurations>
  auto operator()(const Saturated_model& varphi,
                const Point& point,
                R_xlen_t number_types,
                Configurations&&... configurations) const {
    using CountType = std::vector<unsigned long long int>;
    CountType count_vector(number_types);

    for_each_container([&count_vector, &point, &varphi, saturation = varphi.get_saturation(), &configurations...](const auto& current_point) {
      if(count_vector[get_type(current_point)] < saturation && apply_potential(varphi, current_point, point) > 0) {
        ++count_vector[get_type(current_point)];
      }
    }, std::forward<Configurations>(configurations)...);
    // TODO: Ideally, use if constexpr
    if(!Approximate) {
      for_each_container([&count_vector, &point, &varphi, saturation = varphi.get_saturation(), configurations...](const auto& current_point) {
        unsigned long long int count(0);
        for_each_container([&point, &count, &current_point, &varphi, saturation](const auto& other_point) {
          if(!is_equal(other_point, current_point) && get_type(other_point) == get_type(point) && count < saturation && apply_potential(varphi, other_point, current_point) > 0) {
            ++count;
          }
        }, configurations...);
        count_vector[get_type(current_point)] -= count;
        if(!is_equal(point, current_point) && count < saturation && apply_potential(varphi, point, current_point) > 0) {
          ++count;
        }
        count_vector[get_type(current_point)] += count;
      }, std::forward<Configurations>(configurations)...);
    }
    std::vector<double> dispersion(number_types);
    using size_t = typename decltype(dispersion)::size_type;
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      dispersion[i] = static_cast<double>(count_vector[i]);
    }
    return dispersion;
  }
};

template<bool Approximate>
class compute_dispersion<Approximate, dispersionMethod::nonincreasing> {
private:
  using CountType = std::vector<std::vector<double>>;

  template<typename Point, typename Other>
  static void update_count(const Saturated_model& varphi, typename CountType::value_type& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(count.size() < varphi.get_saturation()) {
      count.emplace_back(sq);
      std::push_heap(count.begin(), count.end());
    } else if(sq < count[0]) {
      count.emplace_back(sq);
      std::pop_heap(count.begin(), count.end());
      count.pop_back();
    }
  }
public:
  compute_dispersion() {}

  template<typename Point, typename... Configurations>
  auto operator()(const Saturated_model& varphi,
                const Point& point,
                R_xlen_t number_types,
                Configurations&&... configurations) const {
    CountType count_vector(number_types);

    for_each_container([&count_vector, &point, &varphi](const auto& current_point) {
      compute_dispersion::update_count(varphi, count_vector[get_type(current_point)], current_point, point);
    }, std::forward<Configurations>(configurations)...);

    // TODO: Ideally, use if constexpr
    if(!Approximate) {
      for_each_container([&count_vector, &point, &varphi, configurations...](const auto& current_point) {
        std::vector<double> count{};
        for_each_container([&point, &count, &current_point, &varphi](const auto& other_point) {
          if(!is_equal(other_point, current_point) && get_type(other_point) == get_type(point)) {
            compute_dispersion::update_count(varphi, count, current_point, other_point);
          }
        }, configurations...);
        for(auto c: count) {
          count_vector[get_type(current_point)].emplace_back(-c);
        }
        compute_dispersion::update_count(varphi, count, current_point, point);
        // TODO: Might be able to simplify this and above--only one element of count changes in the computation above.
        for(auto c: count) {
          count_vector[get_type(current_point)].emplace_back(c);
        }
      }, std::forward<Configurations>(configurations)...);
    }
    std::vector<double> dispersion(number_types);
    using size_t = typename decltype(dispersion)::size_type;
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      double d(0.);
      for(auto c: count_vector[i]) {
        // TODO: Don't like this, change?
        if(c < 0) {
          d -= varphi.apply(-c, i, get_type(point));
        } else {
          d += varphi.apply(c, i, get_type(point));
        }
      }
      dispersion[i] = d;
      // dispersion[i] = std::accumulate(count_vector[i].begin(), count_vector[i].end(), 0., [&varphi, i, &point](double count, const auto& val) {
      //   return count + varphi.apply(val, i, get_type(point));
      // });
    }
    return dispersion;
  }
};

template<bool Approximate>
class compute_dispersion<Approximate, dispersionMethod::nonincreasing_after_lower_endpoint> {
private:
  using CountType = std::vector<std::vector<double>>;

  template<typename Point, typename Other>
  static void update_count(const Saturated_model& varphi, typename CountType::value_type& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(sq >= varphi.get_square_lower_endpoint(get_type(point), get_type(other))) {
      if(count.size() < varphi.get_saturation()) {
        count.emplace_back(sq);
        std::push_heap(count.begin(), count.end());
      } else if(sq < count[0]) {
        count.emplace_back(sq);
        std::pop_heap(count.begin(), count.end());
        count.pop_back();
      }
    }
  }
public:
  compute_dispersion() {}

  template<typename Point, typename... Configurations>
  auto operator()(const Saturated_model& varphi,
                const Point& point,
                R_xlen_t number_types,
                Configurations&&... configurations) const {
    CountType count_vector(number_types);

    for_each_container([&count_vector, &point, &varphi](const auto& current_point) {
      compute_dispersion::update_count(varphi, count_vector[get_type(current_point)], current_point, point);
    }, std::forward<Configurations>(configurations)...);

    // TODO: Ideally, use if constexpr
    if(!Approximate) {
      for_each_container([&count_vector, &point, &varphi, configurations...](const auto& current_point) {
        std::vector<double> count{};
        for_each_container([&point, &count, &current_point, &varphi](const auto& other_point) {
          if(!is_equal(other_point, current_point) && get_type(other_point) == get_type(point)) {
            compute_dispersion::update_count(varphi, count, current_point, other_point);
          }
        }, configurations...);
        for(const auto& c: count) {
          count_vector[get_type(current_point)].emplace_back(-c);
        }
        compute_dispersion::update_count(varphi, count, current_point, point);
        // TODO: Might be able to simplify this and above--only one element of count changes in the computation above.
        for(const auto& c: count) {
          count_vector[get_type(current_point)].emplace_back(c);
        }
      }, std::forward<Configurations>(configurations)...);
    }
    std::vector<double> dispersion(number_types);
    using size_t = typename decltype(dispersion)::size_type;
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      double d(0.);
      for(auto c: count_vector[i]) {
        // TODO: Don't like this, change?
        if(c < 0) {
          d -= varphi.apply(-c, i, get_type(point));
        } else {
          d += varphi.apply(c, i, get_type(point));
        }
      }
      dispersion[i] = d;
      // dispersion[i] = std::accumulate(count_vector[i].begin(), count_vector[i].end(), 0., [&varphi, i, &point](double count, const auto& val) {
      //   return count + varphi.apply(val, i, get_type(point));
      // });
    }
    return dispersion;
  }
};

template<bool Approximate>
class compute_dispersion<Approximate, dispersionMethod::generic> {
private:
  using CountType = std::vector<std::vector<double>>;

  template<typename Point, typename Other>
  static void update_count(const Saturated_model& varphi, typename CountType::value_type& count, const Point& point, const Other& other) {
    const auto disp(apply_potential(varphi, other, point));
    if(count.size() < varphi.get_saturation()) {
      count.emplace_back(disp);
      std::push_heap(count.begin(), count.end(), std::greater<double>{});
    } else if(disp > count[0]) {
      count.emplace_back(disp);
      std::pop_heap(count.begin(), count.end(), std::greater<double>{});
      count.pop_back();
    }
  }
public:
  compute_dispersion() {}

  template<typename Point, typename... Configurations>
  auto operator()(const Saturated_model& varphi,
                const Point& point,
                R_xlen_t number_types,
                Configurations&&... configurations) const {
    CountType count_vector(number_types);

    for_each_container([&count_vector, &point, &varphi](const auto& current_point) {
      compute_dispersion::update_count(varphi, count_vector[get_type(current_point)], point, current_point);
    }, std::forward<Configurations>(configurations)...);

    // TODO: Ideally, use if constexpr
    if(!Approximate) {
      for_each_container([&count_vector, &point, &varphi, configurations...](const auto& current_point) {
        std::vector<double> count{};
        for_each_container([&point, &count, &current_point, &varphi](const auto& other_point) {
          if(!is_equal(other_point, current_point) && get_type(other_point) == get_type(point)) {
            compute_dispersion::update_count(varphi, count, other_point, current_point);
          }
        }, configurations...);
        count_vector[get_type(current_point)].emplace_back(-std::accumulate(count.begin(), count.end(), 0.));
        compute_dispersion::update_count(varphi, count, point, current_point);
        count_vector[get_type(current_point)].emplace_back(std::accumulate(count.begin(), count.end(), 0.));
      }, std::forward<Configurations>(configurations)...);
    }

    std::vector<double> dispersion(number_types);
    using size_t = typename decltype(dispersion)::size_type;
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      dispersion[i] = std::accumulate(count_vector[i].begin(), count_vector[i].end(), 0.);
    }
    return dispersion;
  }
};

} // namespace detail

template<bool Approximate = false, typename Point, typename... Configurations>
inline auto compute_dispersion(const Saturated_model& model,
                               const Point& point,
                               R_xlen_t number_types,
                               Configurations&&... configurations) {
  if(model.is_two_valued()) {
    return detail::compute_dispersion<Approximate, detail::dispersionMethod::two_values>{}(model, point, number_types, std::forward<Configurations>(configurations)...);
  } else if(model.is_nonincreasing()) {
    return detail::compute_dispersion<Approximate, detail::dispersionMethod::nonincreasing>{}(model, point, number_types, std::forward<Configurations>(configurations)...);
  } else if(model.is_nonincreasing_after_lower_endpoint()) {
    return detail::compute_dispersion<Approximate, detail::dispersionMethod::nonincreasing_after_lower_endpoint>{}(model, point, number_types, std::forward<Configurations>(configurations)...);
  } else {
    return detail::compute_dispersion<Approximate, detail::dispersionMethod::generic>{}(model, point, number_types, std::forward<Configurations>(configurations)...);
  }
}

// Note: Use inheritance to benefit from EBO.
template<typename Varphi>
class New_saturated_model: public Saturated_model, public Varphi {
public:
  bool is_nonincreasing() const {
    return Varphi::is_nonincreasing;
  }
  bool is_nonincreasing_after_lower_endpoint() const {
    return Varphi::is_nonincreasing_after_lower_endpoint;
  }
  bool is_two_valued() const {
    return Varphi::is_two_valued;
  }

  template<typename... Args>
  New_saturated_model(unsigned long long int saturation, Args&&... args):
    Varphi(std::forward<Args>(args)...),
    saturation_(saturation) {}

  double apply(double normalized_square_distance, int i, int j) const override {
    return Varphi::apply(normalized_square_distance, i, j);
  }
  // TODO: This doesn't return the correct result.
  double get_maximum() const override {
    return static_cast<double>(saturation_);
  }

  unsigned long long int get_saturation() const override {
    return saturation_;
  };

  double get_square_lower_endpoint(int i, int j) const override {
    return Varphi::get_square_lower_endpoint(i, j);
  }
private:
  unsigned long long int saturation_;
};

const constexpr char* const short_range_models[] = {
  "exponential",
  "square_exponential",
  "bump",
  "square_bump",
  "Geyer",
  "linear",
  "Dixon"
};

inline std::unique_ptr<Saturated_model> get_dispersion_from_string(Rcpp::CharacterVector model, Rcpp::NumericMatrix radius, unsigned long long int saturation) {
  const auto model_string(model[0]);
  if(model_string == short_range_models[0]) {
    return std::make_unique<New_saturated_model<potentials::Exponential>>(saturation, radius);
  } else if(model_string == short_range_models[1]) {
    return std::make_unique<New_saturated_model<potentials::Square_exponential>>(saturation, radius);
  } else if(model_string == short_range_models[2]) {
    return std::make_unique<New_saturated_model<potentials::Bump>>(saturation, radius);
  } else if(model_string == short_range_models[3]) {
    return std::make_unique<New_saturated_model<potentials::Square_bump>>(saturation, radius);
  } else if(model_string == short_range_models[4]) {
    return std::make_unique<New_saturated_model<potentials::Strauss>>(saturation, radius);
  } else if(model_string == short_range_models[5]) {
    return std::make_unique<New_saturated_model<potentials::Linear>>(saturation, radius);
  } /*else if(model_string == short_range_models[6]) {
    return f(Dixon_model(saturation));
  } */else {
    Rcpp::stop("Incorrect model entered. A call to show_short_range_models() will show you the available choices.\n");
  }
}

const constexpr char* const medium_range_models[] = {
  "square_exponential",
  "half_square_exponential",
  "Geyer",
  "linear",
  "exponential"
};

inline std::unique_ptr<Saturated_model> get_medium_range_dispersion_from_string(Rcpp::CharacterVector model, Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range, unsigned long long int saturation) {
  const auto model_string(model[0]);
  if(model_string == medium_range_models[0]) {
    return std::make_unique<New_saturated_model<potentials::Medium_range_square_exponential>>(saturation, medium_range, long_range);
  } else if(model_string == medium_range_models[1]) {
    return std::make_unique<New_saturated_model<potentials::Medium_range_half_square_exponential>>(saturation, medium_range, long_range);
  } else if(model_string == medium_range_models[2]) {
    return std::make_unique<New_saturated_model<potentials::Medium_range_Geyer>>(saturation, medium_range, long_range);
  } else if(model_string == medium_range_models[3]) {
    return std::make_unique<New_saturated_model<potentials::Medium_range_linear>>(saturation, medium_range, long_range);
  } else if(model_string == medium_range_models[4]) {
    return std::make_unique<New_saturated_model<potentials::Medium_range_exponential>>(saturation, medium_range, long_range);
  } else {
    Rcpp::stop("Incorrect model entered. A call to show_medium_range_models() will show you the available choices.\n");
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SATURATED_MODEL
