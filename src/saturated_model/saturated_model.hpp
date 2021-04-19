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
#include <memory> // std::shared_ptr
#include <utility> // std::forward, std::move
#include <vector> // std::vector

namespace ppjsdm {

const constexpr char* const short_range_models[] = {
  "exponential",
  "square_exponential",
  "bump",
  "square_bump",
  "Geyer",
  "linear"
};

const constexpr char* const medium_range_models[] = {
  "square_exponential",
  "half_square_exponential",
  "Geyer",
  "linear",
  "half_exponential",
  "exponential",
  "bump",
  "square_bump",
  "tanh"
};

class Saturated_model {
public:
  Saturated_model(Rcpp::CharacterVector model, Rcpp::NumericMatrix radius, unsigned long long int saturation):
  object_(make_short_range_object(model, radius)),
  saturation_(saturation) {}

  Saturated_model(Rcpp::CharacterVector model, Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range, unsigned long long int saturation):
  object_(make_medium_range_object(model, medium_range, long_range)),
  saturation_(saturation) {}

  bool is_nonincreasing() const {
    return object_->is_nonincreasing();
  }

  bool is_nonincreasing_after_lower_endpoint() const {
    return object_->is_nonincreasing_after_lower_endpoint();
  }

  bool is_two_valued() const {
    return object_->is_two_valued();
  }

  double apply(double normalized_square_distance, int i, int j) const {
    return object_->apply(normalized_square_distance, i, j);
  }

  double get_square_lower_endpoint(int i, int j) const {
    return object_->get_square_lower_endpoint(i, j);
  }

  unsigned long long int get_saturation() const {
    return saturation_;
  }

  double get_maximum() const {
    return static_cast<double>(saturation_);
  }

private:
  struct Concept {
    virtual ~Concept() {}
    virtual bool is_nonincreasing() const = 0;
    virtual bool is_nonincreasing_after_lower_endpoint() const = 0;
    virtual bool is_two_valued() const = 0;
    virtual double apply(double normalized_square_distance, int i, int j) const = 0;
    virtual double get_square_lower_endpoint(int i, int j) const = 0;
  };

  template<typename Varphi>
  class Concrete_model: public Concept, private Varphi {
  public:
    bool is_nonincreasing() const override {
      return Varphi::is_nonincreasing;
    }

    bool is_nonincreasing_after_lower_endpoint() const override {
      return Varphi::is_nonincreasing_after_lower_endpoint;
    }

    bool is_two_valued() const override {
      return Varphi::is_two_valued;
    }

    template<typename... Args>
    explicit Concrete_model(Args&&... args):
      Varphi(std::forward<Args>(args)...) {}

    double apply(double normalized_square_distance, int i, int j) const override {
      return Varphi::apply(normalized_square_distance, i, j);
    }

    double get_square_lower_endpoint(int i, int j) const override {
      return Varphi::get_square_lower_endpoint(i, j);
    }
  };

  static std::shared_ptr<const Concept> make_short_range_object(Rcpp::CharacterVector model, Rcpp::NumericMatrix radius) {
    const auto model_string(model[0]);
    if(model_string == short_range_models[0]) {
      return std::make_shared<Concrete_model<potentials::Exponential>>(radius);
    } else if(model_string == short_range_models[1]) {
      return std::make_shared<Concrete_model<potentials::Square_exponential>>(radius);
    } else if(model_string == short_range_models[2]) {
      return std::make_shared<Concrete_model<potentials::Bump>>(radius);
    } else if(model_string == short_range_models[3]) {
      return std::make_shared<Concrete_model<potentials::Square_bump>>(radius);
    } else if(model_string == short_range_models[4]) {
      return std::make_shared<Concrete_model<potentials::Strauss>>(radius);
    } else if(model_string == short_range_models[5]) {
      return std::make_shared<Concrete_model<potentials::Linear>>(radius);
    } else {
      Rcpp::stop("Incorrect model entered. A call to show_short_range_models() will show you the available choices.\n");
    }
  }

  static std::shared_ptr<const Concept> make_medium_range_object(Rcpp::CharacterVector model, Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range) {
    const auto model_string(model[0]);
    if(model_string == medium_range_models[0]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_square_exponential>>(medium_range, long_range);
    } else if(model_string == medium_range_models[1]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_half_square_exponential>>(medium_range, long_range);
    } else if(model_string == medium_range_models[2]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_Geyer>>(medium_range, long_range);
    } else if(model_string == medium_range_models[3]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_linear>>(medium_range, long_range);
    } else if(model_string == medium_range_models[4]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_half_exponential>>(medium_range, long_range);
    } else if(model_string == medium_range_models[5]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_exponential>>(medium_range, long_range);
    } else if(model_string == medium_range_models[6]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_bump>>(medium_range, long_range);
    } else if(model_string == medium_range_models[7]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_square_bump>>(medium_range, long_range);
    } else if(model_string == medium_range_models[8]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_tanh>>(medium_range, long_range);
    } else {
      Rcpp::stop("Incorrect model entered. A call to show_medium_range_models() will show you the available choices.\n");
    }
  }

  std::shared_ptr<const Concept> object_;
  unsigned long long int saturation_;
};

template<typename Point, typename Other>
inline auto apply_potential(const Saturated_model& potential, const Point& point, const Other& other) {
  return potential.apply(normalized_square_distance(point, other), get_type(point), get_type(other));
}

namespace detail {

enum class dispersionMethod {two_values, nonincreasing, nonincreasing_after_lower_endpoint, generic};

template<dispersionMethod Method>
struct compute_dispersion_implementation;

template<>
struct compute_dispersion_implementation<dispersionMethod::two_values> {
  using ValueType = unsigned long long int;

  template<typename Point, typename Vector, typename Integer>
  static auto add_count_to_dispersion(const Saturated_model& varphi,
                                      const Vector& count_vector,
                                      const Point&,
                                      Integer i) {
    const auto count(count_vector[i]);
    return count - static_cast<int>(count > varphi.get_saturation());
  }

  template<typename Point, typename Other>
  static void update_count(const Saturated_model& varphi, ValueType& count, const Point& point, const Other& other) {
    if(count < varphi.get_saturation() + 1 && apply_potential(varphi, point, other) > 0.) {
      ++count;
    }
  }

  template<typename Point, typename Other>
  static void update_count_nonsaturated(const Saturated_model& varphi, ValueType& count, const Point& point, const Other& other) {
    if(apply_potential(varphi, point, other) > 0.) {
      ++count;
    }
  }

  template<bool CountedAlready, typename FloatType, typename Point, typename Other>
  static void add_delta(const Saturated_model& varphi, FloatType& dispersion, const ValueType& count, const Point& point, const Other& other) {
    if(apply_potential(varphi, point, other) > 0.) {
      if(count < varphi.get_saturation() + static_cast<int>(CountedAlready)) {
        ++dispersion;
      }
    }
  }
};

template<>
struct compute_dispersion_implementation<dispersionMethod::nonincreasing> {
  using ValueType = std::vector<double>;

  template<typename Point, typename Vector, typename Integer>
  static auto add_count_to_dispersion(const Saturated_model& varphi,
                                      const Vector& count_vector,
                                      const Point& point,
                                      Integer i) {
    return std::accumulate(count_vector[i].begin() + static_cast<int>(count_vector[i].size() > varphi.get_saturation()),
                           count_vector[i].end(),
                           0., [&varphi, i, &point](double count, const auto& val) {
      return count + varphi.apply(val, i, get_type(point));
    });
  }

  template<typename Point, typename Other>
  static void update_count(const Saturated_model& varphi, ValueType& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(count.size() < varphi.get_saturation() + 1) {
      count.emplace_back(sq);
      std::push_heap(count.begin(), count.end());
    } else if(sq < count[0]) {
      count.emplace_back(sq);
      std::pop_heap(count.begin(), count.end());
      count.pop_back();
    }
  }

  template<typename Point, typename Other>
  static void update_count_nonsaturated(const Saturated_model&, ValueType& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    count.emplace_back(sq);
    std::push_heap(count.begin(), count.end());
  }

  template<bool CountedAlready, typename FloatType, typename Point, typename Other>
  static void add_delta(const Saturated_model& varphi, FloatType& dispersion, const ValueType& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(count.size() < varphi.get_saturation()) {
      dispersion += varphi.apply(sq, get_type(point), get_type(other));
    } else if(count.size() == varphi.get_saturation()) {
      // TODO: if constexpr
      if(CountedAlready) {
        dispersion += varphi.apply(sq, get_type(point), get_type(other));
      } else if(sq < count[0]) {
        dispersion += varphi.apply(sq, get_type(point), get_type(other)) - varphi.apply(count[0], get_type(point), get_type(other));
      }
    } else {
      const auto m(std::max(count[1], count[2]));
      // TODO: if constexpr
      if(CountedAlready) {
        if(sq <= m) {
          dispersion += varphi.apply(sq, get_type(point), get_type(other)) - varphi.apply(count[0], get_type(point), get_type(other));
        }
      } else if(sq < m) {
        dispersion += varphi.apply(sq, get_type(point), get_type(other)) - varphi.apply(m, get_type(point), get_type(other));
      }
    }
  }
};

template<>
struct compute_dispersion_implementation<dispersionMethod::nonincreasing_after_lower_endpoint> {
  using ValueType = std::vector<double>;

  template<typename Point, typename Vector, typename Integer>
  static auto add_count_to_dispersion(const Saturated_model& varphi,
                                      const Vector& count_vector,
                                      const Point& point,
                                      Integer i) {
    return std::accumulate(count_vector[i].begin() + static_cast<int>(count_vector[i].size() > varphi.get_saturation()), count_vector[i].end(), 0., [&varphi, i, &point](double count, const auto& val) {
      return count + varphi.apply(val, i, get_type(point));
    });
  }

  template<typename Point, typename Other>
  static void update_count(const Saturated_model& varphi, ValueType& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(sq >= varphi.get_square_lower_endpoint(get_type(point), get_type(other))) {
      if(count.size() < varphi.get_saturation() + 1) {
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

  template<bool CountedAlready, typename FloatType, typename Point, typename Other>
  static void add_delta(const Saturated_model& varphi, FloatType& dispersion, const ValueType& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(sq >= varphi.get_square_lower_endpoint(get_type(point), get_type(other))) {
      if(count.size() < varphi.get_saturation()) {
        dispersion += varphi.apply(sq, get_type(point), get_type(other));
      } else if(count.size() == varphi.get_saturation()) {
        // TODO: if constexpr
        if(CountedAlready) {
          dispersion += varphi.apply(sq, get_type(point), get_type(other));
        } else if(sq < count[0]) {
          dispersion += varphi.apply(sq, get_type(point), get_type(other)) - varphi.apply(count[0], get_type(point), get_type(other));
        }
      } else {
        const auto m(std::max(count[1], count[2]));
        // TODO: if constexpr
        if(CountedAlready) {
          if(sq <= m) {
            dispersion += varphi.apply(sq, get_type(point), get_type(other)) - varphi.apply(count[0], get_type(point), get_type(other));
          }
        } else if(sq < m) {
            dispersion += varphi.apply(sq, get_type(point), get_type(other)) - varphi.apply(m, get_type(point), get_type(other));
        }
      }
    }
  }
};

template<>
struct compute_dispersion_implementation<dispersionMethod::generic> {
  using ValueType = std::vector<double>;

  template<typename Point, typename Vector, typename Integer>
  static auto add_count_to_dispersion(const Saturated_model& varphi,
                                      const Vector& count_vector,
                                      const Point&,
                                      Integer i) {
    return std::accumulate(count_vector[i].begin() + static_cast<int>(count_vector[i].size() > varphi.get_saturation()),
                           count_vector[i].end(),
                           0.);
  }

  template<typename Point, typename Other>
  static void update_count(const Saturated_model& varphi, ValueType& count, const Point& point, const Other& other) {
    const auto disp(apply_potential(varphi, other, point));
    if(count.size() < varphi.get_saturation() + 1) {
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

  template<bool CountedAlready, typename FloatType, typename Point, typename Other>
  static void add_delta(const Saturated_model& varphi, FloatType& dispersion, const ValueType& count, const Point& point, const Other& other) {
    const auto disp(apply_potential(varphi, other, point));
    if(count.size() < varphi.get_saturation()) {
      dispersion += disp;
    } else if(count.size() == varphi.get_saturation()) {
      // TODO: if constexpr
      if(CountedAlready) {
        dispersion += disp;
      } else if(disp > count[0]) {
        dispersion += disp - count[0];
      }
    } else {
      const auto m(std::min(count[1], count[2]));
      // TODO: if constexpr
      if(CountedAlready) {
        if(disp >= m) {
          dispersion += disp - count[0];
        }
      } else if(disp > m) {
        dispersion += disp - m;
      }
    }
  }
};

template<typename AbstractDispersion, int N = 1, typename Vector, typename CountVector, typename Point>
static void add_count_to_dispersion(const Saturated_model& varphi,
                                    Vector& dispersion,
                                    const CountVector& count_vector,
                                    R_xlen_t number_types,
                                    const Point& point) {
  using size_t = typename Vector::size_type;
  using FloatType = typename Vector::value_type;
  for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
    dispersion[i] += static_cast<FloatType>(N) * static_cast<FloatType>(AbstractDispersion::add_count_to_dispersion(varphi, count_vector, point, i));
  }
}

template<typename AbstractDispersion, typename Point, typename... Configurations>
inline auto generic_dispersion_computation(const Saturated_model& varphi,
                                           const Point& point,
                                           R_xlen_t number_types,
                                           Configurations&... configurations) {
  // TODO: Something important to think about:
  // When running gibbsm or vcov.gibbsm, we need the saturation of every point in the configuration w.r.t.
  // each other point in the configuration.
  // In the current way I'm running this, we end up computing saturation of each point multiple times.
  // Large speed gains are possible by finding a better strategy.

  // TODO: The max size of the heap I'm using is always the saturation; reserve.
  using ValueType = typename AbstractDispersion::ValueType;
  using CountType = std::vector<ValueType>;
  using DispersionType = std::vector<double>;

  CountType count_vector(number_types);
  DispersionType dispersion(number_types);
  if(static_cast<decltype(size(configurations...))>(varphi.get_saturation()) >= size(configurations...)) {
    for_each_container([&count_vector, &point, &varphi](const auto& current_point) {
      if(!is_equal(current_point, point)) {
        AbstractDispersion::update_count_nonsaturated(varphi, count_vector[get_type(current_point)], current_point, point);
      }
    }, configurations...);
    add_count_to_dispersion<AbstractDispersion, 2>(varphi, dispersion, count_vector, number_types, point);
  } else {
    for_each_container([&dispersion, &count_vector, &point, &varphi,
                       saturation = varphi.get_saturation(), &configurations...](const auto& current_point) {
                         if(!is_equal(current_point, point)) {
                           ValueType count{};
                           // TODO: In the Geyer case, there's potential to use conditional_for_each_container and break if count == saturation.
                           // However, note that this would make the function different from other variations.
                           for_each_container([&point, &count, &current_point, &varphi, saturation](const auto& other_point) {
                             if(!is_equal(other_point, point) && !is_equal(other_point, current_point) && get_type(other_point) == get_type(point)) {
                               AbstractDispersion::update_count(varphi, count, current_point, other_point);
                             }
                           }, configurations...);
                           AbstractDispersion::update_count(varphi, count_vector[get_type(current_point)], current_point, point);
                           AbstractDispersion::template add_delta<false>(varphi, dispersion[get_type(current_point)], count, current_point, point);
                         }
                       }, configurations...);
    add_count_to_dispersion<AbstractDispersion, 1>(varphi, dispersion, count_vector, number_types, point);
  }
  return dispersion;
}

template<typename AbstractDispersion, typename Configuration, typename OtherConfiguration>
inline auto generic_dispersion_computation(const Saturated_model& varphi,
                                           R_xlen_t number_types,
                                           const Configuration& configuration,
                                           const OtherConfiguration& other_configuration) {
  using ValueType = typename AbstractDispersion::ValueType;
  using CountType = std::vector<ValueType>;
  using DispersionType = std::vector<double>;

  const auto configuration_size(size(configuration));
  using size_t = std::remove_cv_t<decltype(size(configuration))>;

  std::vector<CountType> count_vector(configuration_size);
  for(size_t i(0); i < configuration_size; ++i) {
    count_vector[i] = CountType(number_types);
    for(size_t j(0); j < configuration_size; ++j) {
      if(i != j) {
        AbstractDispersion::update_count(varphi, count_vector[i][get_type(configuration[j])], configuration[i], configuration[j]);
      }
    }
  }

  const auto other_configuration_size(size(other_configuration));
  std::vector<DispersionType> dispersion(configuration_size + other_configuration_size);

  for(size_t i(0); i < configuration_size; ++i) {
    dispersion[i] = DispersionType(number_types);
    for(size_t j(0); j < configuration_size; ++j) {
      if(i != j) {
        AbstractDispersion::template add_delta<true>(varphi, dispersion[i][get_type(configuration[j])], count_vector[j][get_type(configuration[i])], configuration[j], configuration[i]);
      }
    }
    add_count_to_dispersion<AbstractDispersion, 1>(varphi, dispersion[i], count_vector[i], number_types, configuration[i]);
  }

  for(size_t i(0); i < other_configuration_size; ++i) {
    dispersion[configuration_size + i] = DispersionType(number_types);
    CountType count_point(number_types);
    for(size_t j(0); j < configuration_size; ++j) {
      AbstractDispersion::template add_delta<false>(varphi, dispersion[configuration_size + i][get_type(configuration[j])], count_vector[j][get_type(other_configuration[i])], configuration[j], other_configuration[i]);
      AbstractDispersion::update_count(varphi, count_point[get_type(configuration[j])], other_configuration[i], configuration[j]);
    }
    add_count_to_dispersion<AbstractDispersion, 1>(varphi, dispersion[configuration_size + i], count_point, number_types, other_configuration[i]);
  }

  return dispersion;
}

} // namespace detail

template<typename... Configurations>
inline auto compute_dispersion_index(const Saturated_model& model,
                                     R_xlen_t number_types,
                                     Configurations&&... configurations) {
  if(model.is_two_valued()) {
    return detail::generic_dispersion_computation<detail::compute_dispersion_implementation<detail::dispersionMethod::two_values>>(model, number_types, std::forward<Configurations>(configurations)...);
  } else if(model.is_nonincreasing()) {
    return detail::generic_dispersion_computation<detail::compute_dispersion_implementation<detail::dispersionMethod::nonincreasing>>(model, number_types, std::forward<Configurations>(configurations)...);
  } else if(model.is_nonincreasing_after_lower_endpoint()) {
    return detail::generic_dispersion_computation<detail::compute_dispersion_implementation<detail::dispersionMethod::nonincreasing_after_lower_endpoint>>(model, number_types, std::forward<Configurations>(configurations)...);
  } else {
    return detail::generic_dispersion_computation<detail::compute_dispersion_implementation<detail::dispersionMethod::generic>>(model, number_types, std::forward<Configurations>(configurations)...);
  }
}

template<typename Point, typename... Configurations>
inline auto compute_dispersion(const Saturated_model& model,
                               const Point& point,
                               R_xlen_t number_types,
                               Configurations&&... configurations) {
  if(model.is_two_valued()) {
    return detail::generic_dispersion_computation<detail::compute_dispersion_implementation<detail::dispersionMethod::two_values>>(model, point, number_types, std::forward<Configurations>(configurations)...);
  } else if(model.is_nonincreasing()) {
    return detail::generic_dispersion_computation<detail::compute_dispersion_implementation<detail::dispersionMethod::nonincreasing>>(model, point, number_types, std::forward<Configurations>(configurations)...);
  } else if(model.is_nonincreasing_after_lower_endpoint()) {
    return detail::generic_dispersion_computation<detail::compute_dispersion_implementation<detail::dispersionMethod::nonincreasing_after_lower_endpoint>>(model, point, number_types, std::forward<Configurations>(configurations)...);
  } else {
    return detail::generic_dispersion_computation<detail::compute_dispersion_implementation<detail::dispersionMethod::generic>>(model, point, number_types, std::forward<Configurations>(configurations)...);
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SATURATED_MODEL
