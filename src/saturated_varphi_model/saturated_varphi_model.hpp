#ifndef INCLUDE_PPJSDM_SATURATED_VARPHI
#define INCLUDE_PPJSDM_SATURATED_VARPHI

#include <Rcpp.h>

#include "varphi.hpp"
#include "../configuration/configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"
#include "../point/square_distance.hpp"
#include "../utility/for_each_container.hpp"

#include <algorithm> // std::upper_bound
#include <list> // std::list
#include <utility> // std::forward, std::move
#include <vector> // std::vector

namespace ppjsdm {

template<typename Point, typename Dispersion, typename Select = void>
class Compute_dispersion;

template<typename Point, typename Dispersion>
class Compute_dispersion<Point, Dispersion, std::enable_if_t<has_nonzero_value_v<Dispersion>>> {
private:
  using CountType = std::vector<unsigned long long int>;
public:
  Compute_dispersion(const Point& point,
                     R_xlen_t number_types,
                     const Dispersion& dispersion):
  point_(point),
  dispersion_(dispersion),
  count_positive_types_(number_types) {}

  template<typename Configuration>
  void add_configuration(const Configuration& configuration) {
    update_count(configuration);
  }

  auto compute() const {
    return compute_dispersion(count_positive_types_);
  }

  auto compute_from_state(const CountType& count) const {
    return compute_dispersion(count);
  }

  auto get_state() const {
    return count_positive_types_;
  }
private:
  const Point& point_;
  const Dispersion& dispersion_;
  CountType count_positive_types_;

  template<typename Configuration>
  void update_count(const Configuration& configuration) {
    using size_t = size_t<Configuration>;
    for(size_t i(0); i < size(configuration); ++i) {
      const auto& current_point(configuration[i]);
      if(count_positive_types_[get_type(current_point)] < dispersion_.saturation_ && dispersion_.apply(current_point, point_) > 0) {
        ++count_positive_types_[get_type(current_point)];
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
class Compute_dispersion<Point, Dispersion, std::enable_if_t<!has_nonzero_value_v<Dispersion>>> {
private:
  using CountType = std::vector<std::list<double>>;
public:
  Compute_dispersion(const Point& point,
                     R_xlen_t number_types,
                     const Dispersion& dispersion):
  point_(point),
  dispersion_(dispersion),
  square_distances_(number_types) {}

  template<typename Configuration>
  void add_configuration(const Configuration& configuration) {
    update_count(configuration);
  }

  auto compute() const {
    return compute_dispersion(square_distances_);
  }

  auto compute_from_state(const CountType& count) const {
    return compute_dispersion(count);
  }

  auto get_state() const {
    return square_distances_;
  }
private:
  const Point& point_;
  const Dispersion& dispersion_;
  CountType square_distances_;

  template<typename Configuration>
  void update_count(const Configuration& configuration) {
    using size_t = size_t<Configuration>;
    for(size_t i(0); i < size(configuration); ++i) {
      const auto& current_point(configuration[i]);
      const auto sq(square_distance(current_point, point_));
      auto& current(square_distances_[get_type(current_point)]);
      auto iterator(std::upper_bound(current.begin(), current.end(), sq));
      if(current.size() < dispersion_.saturation_) {
        current.insert(iterator, sq);
      } else if(iterator != current.end()) {
        current.insert(iterator, sq);
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
      for(const auto sq: square_distances_[i]) {
        d += dispersion_.apply(sq, i, point_type);
      }
      dispersion[i] = d;
    }
    return dispersion;
  }
};

// Note: Use inheritance to benefit from EBO.
// TODO: Try to get rid of public inheritance if possible.
template<typename Varphi>
class Saturated_varphi_model: public Varphi {
private:
  template<typename, typename, typename>
  friend class Compute_dispersion;
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
