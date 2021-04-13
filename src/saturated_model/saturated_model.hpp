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
  using CountType = std::vector<unsigned long long int>;

  template<int N = 1, typename Vector, typename Point>
  static void update_dispersion(Vector& dispersion,
                                const CountType& count_vector,
                                R_xlen_t number_types,
                                const Saturated_model&,
                                const Point&) {
    using size_t = typename Vector::size_type;
    using FloatType = typename Vector::value_type;
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      dispersion[i] += static_cast<FloatType>(N) * static_cast<FloatType>(count_vector[i]);
    }
  }

  template<typename Point, typename Other>
  static void update_count(const Saturated_model& varphi, typename CountType::value_type& count, const Point& point, const Other& other) {
    if(count < varphi.get_saturation() && apply_potential(varphi, point, other) > 0.) {
      ++count;
    }
  }

  template<typename Point, typename Other>
  static void update_count_nonsaturated(const Saturated_model& varphi, typename CountType::value_type& count, const Point& point, const Other& other) {
    if(apply_potential(varphi, point, other) > 0.) {
      ++count;
    }
  }

  template<typename Vector, typename Point, typename Other>
  static void add_delta(Vector& dispersion, const Saturated_model& varphi, const typename CountType::value_type& count, const Point& point, const Other& other) {
    if(count < varphi.get_saturation() && apply_potential(varphi, point, other) > 0.) {
      ++dispersion[get_type(point)];
    }
  }
};

template<>
struct compute_dispersion_implementation<dispersionMethod::nonincreasing> {
  using CountType = std::vector<std::vector<double>>;

  template<int N = 1, typename Vector, typename Point>
  static void update_dispersion(Vector& dispersion,
                                const CountType& count_vector,
                                R_xlen_t number_types,
                                const Saturated_model& varphi,
                                const Point& point) {
    using size_t = typename Vector::size_type;
    using FloatType = typename Vector::value_type;
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      dispersion[i] += static_cast<FloatType>(N) * std::accumulate(count_vector[i].begin(), count_vector[i].end(), 0., [&varphi, i, &point](double count, const auto& val) {
        return count + varphi.apply(val, i, get_type(point));
      });
    }
  }

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

  template<typename Point, typename Other>
  static void update_count_nonsaturated(const Saturated_model&, typename CountType::value_type& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    count.emplace_back(sq);
    std::push_heap(count.begin(), count.end());
  }

  template<typename Vector, typename Point, typename Other>
  static void add_delta(Vector& dispersion, const Saturated_model& varphi, const typename CountType::value_type& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(count.size() < varphi.get_saturation()) {
      dispersion[get_type(point)] += varphi.apply(sq, get_type(point), get_type(other));
    } else if(sq < count[0]) {
      dispersion[get_type(point)] += varphi.apply(sq, get_type(point), get_type(other)) - varphi.apply(count[0], get_type(point), get_type(other));
    }
  }
};

template<>
struct compute_dispersion_implementation<dispersionMethod::nonincreasing_after_lower_endpoint> {
  using CountType = std::vector<std::vector<double>>;

  template<int N = 1, typename Vector, typename Point>
  static void update_dispersion(Vector& dispersion,
                                const CountType& count_vector,
                                R_xlen_t number_types,
                                const Saturated_model& varphi,
                                const Point& point) {
    using size_t = typename Vector::size_type;
    using FloatType = typename Vector::value_type;
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      dispersion[i] += static_cast<FloatType>(N) * std::accumulate(count_vector[i].begin(), count_vector[i].end(), 0., [&varphi, i, &point](double count, const auto& val) {
        return count + varphi.apply(val, i, get_type(point));
      });
    }
  }

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

  template<typename Point, typename Other>
  static void update_count_nonsaturated(const Saturated_model& varphi, typename CountType::value_type& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(sq >= varphi.get_square_lower_endpoint(get_type(point), get_type(other))) {
      count.emplace_back(sq);
      std::push_heap(count.begin(), count.end());
    }
  }

  template<typename Vector, typename Point, typename Other>
  static void add_delta(Vector& dispersion, const Saturated_model& varphi, const typename CountType::value_type& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(sq >= varphi.get_square_lower_endpoint(get_type(point), get_type(other))) {
      if(count.size() < varphi.get_saturation()) {
        dispersion[get_type(point)] += varphi.apply(sq, get_type(point), get_type(other));
      } else if(sq < count[0]) {
        dispersion[get_type(point)] += varphi.apply(sq, get_type(point), get_type(other)) - varphi.apply(count[0], get_type(point), get_type(other));
      }
    }
  }
};

template<>
struct compute_dispersion_implementation<dispersionMethod::generic> {
  using CountType = std::vector<std::vector<double>>;

  template<int N = 1, typename Vector, typename Point>
  static void update_dispersion(Vector& dispersion,
                                const CountType& count_vector,
                                R_xlen_t number_types,
                                const Saturated_model&,
                                const Point&) {
    using size_t = typename Vector::size_type;
    using FloatType = typename Vector::value_type;
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      dispersion[i] += static_cast<FloatType>(N) * std::accumulate(count_vector[i].begin(), count_vector[i].end(), 0.);
    }
  }

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

  template<typename Point, typename Other>
  static void update_count_nonsaturated(const Saturated_model& varphi, typename CountType::value_type& count, const Point& point, const Other& other) {
    const auto disp(apply_potential(varphi, other, point));
    count.emplace_back(disp);
    std::push_heap(count.begin(), count.end(), std::greater<double>{});
  }

  template<typename Vector, typename Point, typename Other>
  static void add_delta(Vector& dispersion, const Saturated_model& varphi, const typename CountType::value_type& count, const Point& point, const Other& other) {
    const auto disp(apply_potential(varphi, other, point));
    if(count.size() < varphi.get_saturation()) {
      dispersion[get_type(point)] += disp;
    } else if(disp > count[0]) {
      dispersion[get_type(point)] += disp - count[0];
    }
  }
};

template<typename AbstractDispersion, typename Point, typename... Configurations>
inline auto generic_dispersion_computation(const Saturated_model& varphi,
                                           const Point& point,
                                           R_xlen_t number_types,
                                           Configurations&... configurations) {
  using DispersionType = std::vector<double>;

  typename AbstractDispersion::CountType count_vector(number_types);
  DispersionType dispersion(number_types);
  if(static_cast<decltype(size(configurations...))>(varphi.get_saturation()) >= size(configurations...)) {
    for_each_container([&count_vector, &point, &varphi](const auto& current_point) {
      if(!is_equal(current_point, point)) {
        AbstractDispersion::update_count_nonsaturated(varphi, count_vector[get_type(current_point)], current_point, point);
      }
    }, configurations...);
    AbstractDispersion::template update_dispersion<2>(dispersion, count_vector, number_types, varphi, point);
  } else {
    for_each_container([&dispersion, &count_vector, &point, &varphi,
                       saturation = varphi.get_saturation(), &configurations...](const auto& current_point) {
                         if(!is_equal(current_point, point)) {
                           typename AbstractDispersion::CountType::value_type count{};
                           // TODO: In the Geyer case, there's potential to use conditional_for_each_container and break if count == saturation.
                           // However, note that this would make the function different from other variations.
                           for_each_container([&point, &count, &current_point, &varphi, saturation](const auto& other_point) {
                             if(!is_equal(other_point, point) && !is_equal(other_point, current_point) && get_type(other_point) == get_type(point)) {
                               AbstractDispersion::update_count(varphi, count, current_point, other_point);
                             }
                           }, configurations...);
                           AbstractDispersion::update_count(varphi, count_vector[get_type(current_point)], current_point, point);
                           AbstractDispersion::add_delta(dispersion, varphi, count, current_point, point);
                         }
                       }, configurations...);
    AbstractDispersion::template update_dispersion<1>(dispersion, count_vector, number_types, varphi, point);
  }
  return dispersion;
}

} // namespace detail

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
