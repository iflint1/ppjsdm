#ifndef INCLUDE_PPJSDM_SATURATED_VARPHI
#define INCLUDE_PPJSDM_SATURATED_VARPHI

#include <Rcpp.h>

#include "potential.hpp"
#include "../configuration/configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"
#include "../point/square_distance.hpp"
#include "../utility/for_each_container.hpp"

#include <algorithm> // std::accumulate, std::pop_heap, std::push_heap
#include <functional> // std::greater
#include <utility> // std::forward, std::move
#include <vector> // std::vector

namespace ppjsdm {
namespace detail {

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
    using size_t = size_t<Configuration>;
    for(size_t i(0); i < size(configuration); ++i) {
      const auto& current_point(configuration[i]);
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
class Compute_dispersion_implementation<Point, Dispersion, std::enable_if_t<!has_nonzero_value_v<Dispersion> && Dispersion::is_nonincreasing>> {
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
    using size_t = size_t<Configuration>;
    for(size_t i(0); i < size(configuration); ++i) {
      const auto& current_point(configuration[i]);
      const auto sq(square_distance(current_point, point_));
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
class Compute_dispersion_implementation<Point, Dispersion, std::enable_if_t<!has_nonzero_value_v<Dispersion> && !Dispersion::is_nonincreasing>> {
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
    using size_t = size_t<Configuration>;
    for(size_t i(0); i < size(configuration); ++i) {
      const auto& current_point(configuration[i]);
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

  auto get_state() const {
    return ImplementationBase::count_vector_;
  }
};

// Note: Use inheritance to benefit from EBO.
template<typename Varphi>
class Saturated_varphi_model: public Varphi {
private:
  template<typename, typename, typename>
  friend class detail::Compute_dispersion_implementation;
public:
  template<typename... Args>
  Saturated_varphi_model(unsigned long long int saturation, Args&&... args):
    Varphi(std::forward<Args>(args)...), saturation_(saturation) {}

  template<typename Point>
  auto get_compute_dispersion_object(const Point& point,
                                     R_xlen_t number_types) const {
    return Compute_dispersion<Point, Saturated_varphi_model<Varphi>>(point, number_types, *this);
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

const constexpr char* const models[] = {
  "exponential",
  "square_exponential",
  "bump",
  "square_bump",
  "Geyer"
};

template<typename F>
inline auto call_on_dispersion_model(Rcpp::CharacterVector model, Rcpp::NumericMatrix radius, unsigned long long int saturation, const F& f) {
  const auto model_string(model[0]);
  if(model_string == models[0]) {
    return f(Saturated_varphi_model<varphi::Exponential>(saturation, radius));
  } else if(model_string == models[1]) {
    return f(Saturated_varphi_model<varphi::Square_exponential>(saturation, radius));
  } else if(model_string == models[2]) {
    return f(Saturated_varphi_model<varphi::Bump>(saturation, radius));
  } else if(model_string == models[3]) {
    return f(Saturated_varphi_model<varphi::Square_bump>(saturation, radius));
  } else if(model_string == models[4]) {
    return f(Saturated_varphi_model<varphi::Strauss>(saturation, radius));
  } else {
    Rcpp::stop("Incorrect model entered. A call to show_models() will show you the available choices.\n");
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SATURATED_VARPHI
