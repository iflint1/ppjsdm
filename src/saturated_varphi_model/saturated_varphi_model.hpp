#ifndef INCLUDE_PPJSDM_SATURATED_VARPHI
#define INCLUDE_PPJSDM_SATURATED_VARPHI

#include <Rcpp.h>

#include "varphi.hpp"
#include "../point/point_manipulation.hpp"
#include "../point/square_distance.hpp"
#include "../utility/for_each_container.hpp"

#include <algorithm> // std::upper_bound
#include <list> // std::list
#include <utility> // std::forward, std::move
#include <vector> // std::vector

namespace ppjsdm {

// Note: Use public inheritance to benefit from EBO.
template<typename Varphi>
class Saturated_varphi_model: public Varphi {
public:
  template<typename... Args>
  Saturated_varphi_model(unsigned long long int saturation, Args&&... args):
    Varphi(std::forward<Args>(args)...), saturation_(saturation) {}

  template<typename Point, typename... Configurations>
  auto compute(const Point& point,
               R_xlen_t number_types,
               Configurations&&... configurations) const {
    // TODO: Ideally I'd like to use if constexpr
    // TODO: If I can use if constexpr, I can remove get_nonzero_value below.
    if(has_nonzero_value_v<Varphi>) { // Use faster algorithm in this case
      std::vector<unsigned long long int> count_positive_types(number_types);
      for_each_container([&count_positive_types, this, &point](const auto& current_point) {
        if(count_positive_types[get_type(current_point)] < saturation_ && this->Varphi::apply(current_point, point) > 0) {
          count_positive_types[get_type(current_point)] += 1;
        }
      }, std::forward<Configurations>(configurations)...);

      std::vector<double> dispersion(number_types);
      using size_t = typename decltype(dispersion)::size_type;
      constexpr double nonzero_value(get_nonzero_value<Varphi>());
      for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
        dispersion[i] = nonzero_value * static_cast<double>(count_positive_types[i]);
      }
      return dispersion;
    } else {
      std::vector<std::list<double>> square_distances(number_types);

      // Fill with `saturation_` smallest square distances
      for_each_container([&square_distances, &point, saturation = saturation_](const auto& current_point) {
        const auto sq(square_distance(current_point, point));
        auto& current(square_distances[get_type(current_point)]);
        auto iterator(std::upper_bound(current.begin(), current.end(), sq));
        if(current.size() < saturation) {
          current.insert(iterator, sq);
        } else if(iterator != current.end()) {
          current.insert(iterator, sq);
          current.pop_back();
        }
      }, std::forward<Configurations>(configurations)...);

      // Compute dispersion
      std::vector<double> dispersion(number_types);
      using size_t = typename decltype(dispersion)::size_type;
      const auto point_type(get_type(point));
      for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
        double d(0);
        for(const auto sq: square_distances[i]) {
          d += Varphi::apply(sq, i, point_type);
        }
        dispersion[i] = d;
      }
      return dispersion;
    }
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
