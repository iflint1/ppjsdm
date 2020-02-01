#ifndef INCLUDE_PPJSDM_SATURATED_VARPHI
#define INCLUDE_PPJSDM_SATURATED_VARPHI

#include <Rcpp.h>

#include "varphi.hpp"
#include "../configuration/configuration_manipulation.hpp"
#include "../configuration/get_number_points.hpp"
#include "../point/point_manipulation.hpp"
#include "../utility/im_wrapper.hpp"
#include "../utility/size_t.hpp"

#include <algorithm> // std::max, std::sort
#include <cmath> // std::exp, std::log, std::fabs
#include <type_traits> // std::remove_const, std::remove_reference, std::is_same, std::enable_if
#include <utility> // std::forward
#include <vector> // std::vector

namespace ppjsdm {

// Note: Use public inheritance to benefit from EBO.
template<typename Varphi>
class Saturated_varphi_model_papangelou: public Varphi {
public:
  template<typename... Args>
  Saturated_varphi_model_papangelou(unsigned long long int saturation, Args&&... args):
    Varphi(std::forward<Args>(args)...), saturation_(saturation) {}

  template<typename Configuration, typename Point>
  Rcpp::NumericVector compute(const Configuration& configuration,
                              const Point& point,
                              R_xlen_t number_types,
                              size_t<Configuration> number_points) const {
    std::vector<std::vector<double>> square_distances(number_types);

    // Preallocate enough to account for worst case.
    for(R_xlen_t i(0); i < number_types; ++i) {
      square_distances[i].reserve(number_points);
    }

    // Fill with square distances
    using size_t = size_t<Configuration>;
    const auto point_x(get_x(point));
    const auto point_y(get_y(point));
    for(size_t i(0); i < number_points; ++i) {
      const auto point_i(configuration[i]);
      const auto type_i(get_type(point_i));
      const auto delta_x(get_x(point_i) - point_x);
      const auto delta_y(get_y(point_i) - point_y);
      square_distances[type_i].emplace_back(delta_x * delta_x + delta_y * delta_y);
    }

    // Compute dispersion
    Rcpp::NumericVector dispersion(Rcpp::no_init(number_types));
    const auto point_type(get_type(point));
    for(R_xlen_t i(0); i < number_types; ++i) {
      auto current_square_distances(square_distances[i]);
      std::sort(current_square_distances.begin(), current_square_distances.end());
      const auto points_to_consider(current_square_distances.size() < saturation_
                                      ? current_square_distances.size()
                                      : saturation_);
      double disp(0);
      for(decltype(saturation_) j(0); j < points_to_consider; ++j) {
        disp += Varphi::apply(current_square_distances[j], i, point_type);
      }
      if(points_to_consider > 0) {
        disp /= static_cast<double>(points_to_consider);
      }
      dispersion[i] = disp;
    }
    return dispersion;
  }

  template<typename Window>
  double get_maximum(const Window& window) const {
    return Varphi::get_maximum(window);
  }
private:
  unsigned long long int saturation_;
};

// Note: Use public inheritance to benefit from EBO.
template<typename Varphi>
class Mean_varphi_model_papangelou: public Varphi {
public:
  template<typename... Args>
  Mean_varphi_model_papangelou(Args&&... args): Varphi(std::forward<Args>(args)...) {}

  template<typename Configuration, typename Point>
  Rcpp::NumericVector compute(const Configuration& configuration,
                              const Point& point,
                              R_xlen_t number_types,
                              size_t<Configuration> number_points) const {
    using size_t = size_t<Configuration>;

    // dispersion and count_types are automatically 0-initialized
    Rcpp::NumericVector dispersion(number_types);
    std::vector<size_t> count_types(number_types);

    const auto type_point(get_type(point));
    const auto x_point(get_x(point));
    const auto y_point(get_y(point));
    for(size_t i(0); i < number_points; ++i) {
      const auto point_i(configuration[i]);
      const auto type_i(get_type(point_i));
      ++count_types[type_i];
      const auto delta_x(get_x(point_i) - x_point);
      const auto delta_y(get_y(point_i) - y_point);
      dispersion[type_i] += Varphi::apply(delta_x * delta_x + delta_y * delta_y, type_i, type_point);
    }

    for(R_xlen_t i(0); i < number_types; ++i) {
      const auto count_types_i(count_types[i]);
      if(count_types_i > 0) {
        dispersion[i] /= static_cast<double>(count_types_i);
      }
    }

    return dispersion;
  }

  template<typename Window>
  double get_maximum(const Window& window) const {
    return Varphi::get_maximum(window);
  }
};

template<typename Dispersion, typename Lambda, typename Alpha, typename Coefs>
class Exponential_family_model: public Dispersion {
public:
  template<typename D, std::enable_if_t<std::is_same<Dispersion, std::remove_reference_t<D>>::value>* = nullptr>
  Exponential_family_model(const Lambda& lambda,
                           const Alpha& alpha,
                           const Coefs& coefs,
                           Rcpp::List covariates,
                           D&& dispersion):
  Dispersion(std::forward<D>(dispersion)),
    lambda_(lambda),
    alpha_(alpha),
    coefs_(coefs),
    covariates_(covariates) {}

  template<typename... Args>
  Exponential_family_model(const Lambda& lambda,
                           const Alpha& alpha,
                           const Coefs& coefs,
                           Rcpp::List covariates,
                           Args&&... args):
    Exponential_family_model(lambda, alpha, coefs, covariates, Dispersion(std::forward<Args>(args)...)) {}

  template<typename Configuration, typename Point>
  double compute_papangelou(const Configuration& configuration,
                            const Point& point,
                            R_xlen_t number_types) const {
    const auto dispersion(Dispersion::compute(configuration, point, number_types, size(configuration)));

    double inner_product(0);
    const auto point_type(get_type(point));
    for(R_xlen_t i(0); i < number_types; ++i) {
      inner_product += alpha_(i, point_type) * dispersion[i];
    }
    for(decltype(covariates_.size()) i(0); i < covariates_.size(); ++i) {
      inner_product += coefs_(i, point_type) * covariates_[i](get_x(point), get_y(point));
    }
    return lambda_[point_type] * std::exp(inner_product);
  }

  template<typename Configuration, typename Point, typename Other>
  double compute_papangelou_conditional_on_value(const Configuration& configuration,
                            const Point& point,
                            const Other& other,
                            double,
                            R_xlen_t number_types) const {
    Configuration copy(configuration);
    for(const auto& p: other) {
      add_point(copy, p);
    }
    return compute_papangelou(copy, point, number_types);
  }

  template<typename Window>
  auto get_dominating_intensity(const Window& window, R_xlen_t number_types) const {
    std::vector<double> dominating_intensity(number_types);
    for(R_xlen_t i(0); i < number_types; ++i) {
      double inner_product(0);
      for(R_xlen_t j(0); j < number_types; ++j) {
        inner_product += std::fabs(alpha_(j, i)) * Dispersion::get_maximum(window);
      }
      for(decltype(covariates_.size()) j(0); j < covariates_.size(); ++j) {
        const auto bounds(covariates_[j].bounds());
        inner_product += std::fabs(coefs_(j, i)) * std::max(std::fabs(bounds.first), std::fabs(bounds.second));
      }
      const auto value(lambda_[i] * std::exp(inner_product));
      if(std::isinf(value)) {
        Rcpp::stop("Infinite value obtained as the bound to the Papangelou intensity.");
      }
      dominating_intensity[i] = value;
    }
    return dominating_intensity;
  }
private:
  Lambda lambda_;
  Alpha alpha_;
  Coefs coefs_;
  Im_list_wrapper covariates_;
};

const constexpr char* const models[] = {
  "exponential",
  "square_exponential",
  "Strauss",
  "saturated_exponential",
  "saturated_square_exponential",
  "Geyer"
};

template<typename F>
inline auto call_on_papangelou(Rcpp::CharacterVector model, Rcpp::NumericMatrix radius, unsigned long long int saturation, const F& f) {
  const auto model_string(model[0]);
  if(model_string == models[0]) {
    return f(Mean_varphi_model_papangelou<varphi::Exponential>(radius));
  } else if(model_string == models[1]) {
    return f(Mean_varphi_model_papangelou<varphi::Square_exponential>(radius));
  } else if(model_string == models[2]) {
    return f(Mean_varphi_model_papangelou<varphi::Strauss>(radius));
  } else if(model_string == models[3]) {
    return f(Saturated_varphi_model_papangelou<varphi::Exponential>(saturation, radius));
  } else if(model_string == models[4]) {
    return f(Saturated_varphi_model_papangelou<varphi::Square_exponential>(saturation, radius));
  } else if(model_string == models[5]) {
    return f(Saturated_varphi_model_papangelou<varphi::Strauss>(saturation, radius));
  } else {
    Rcpp::stop("Incorrect model entered. A call to show_models() will show you the available choices.\n");
  }
}

template<typename F, typename Lambda>
inline auto call_on_model(Rcpp::CharacterVector model,
                          Rcpp::NumericMatrix alpha,
                          const Lambda& lambda,
                          Rcpp::NumericMatrix coefs,
                          Rcpp::List covariates,
                          Rcpp::NumericMatrix radius,
                          unsigned long long int saturation,
                          const F& f) {
  return call_on_papangelou(model, radius, saturation, [&alpha, &lambda, &f, &coefs, covariates](auto&& varphi) {
    using Model_type = Exponential_family_model<std::remove_cv_t<std::remove_reference_t<decltype(varphi)>>,
                                       Lambda,
                                       Rcpp::NumericMatrix,
                                       Rcpp::NumericMatrix>;
    const Model_type model(lambda, alpha, coefs, covariates, std::forward<decltype(varphi)>(varphi));
    return f(model);
  });
}


} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SATURATED_VARPHI
