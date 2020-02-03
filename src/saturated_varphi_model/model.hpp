#ifndef INCLUDE_PPJSDM_SATURATED_VARPHI
#define INCLUDE_PPJSDM_SATURATED_VARPHI

#include <Rcpp.h>

#include "varphi.hpp"
#include "../configuration/configuration_manipulation.hpp"
#include "../configuration/get_number_points.hpp"
#include "../point/point_manipulation.hpp"
#include "../point/square_distance.hpp"
#include "../utility/for_each_container.hpp"
#include "../utility/im_wrapper.hpp"
#include "../utility/size_t.hpp"

#include <algorithm> // std::max, std::min, std::upper_bound
#include <cmath> // std::exp, std::log, std::fabs, std::isinf
#include <list> // std::list
#include <type_traits> // std::remove_const, std::remove_reference, std::is_same, std::enable_if
#include <utility> // std::forward, std::move
#include <vector> // std::vector

namespace ppjsdm {

// Note: Use public inheritance to benefit from EBO.
template<typename Varphi>
class Saturated_varphi_model_papangelou: public Varphi {
public:
  template<typename... Args>
  Saturated_varphi_model_papangelou(unsigned long long int saturation, Args&&... args):
    Varphi(std::forward<Args>(args)...), saturation_(saturation) {}

  template<typename... Configurations, typename Point>
  Rcpp::NumericVector compute(const Point& point,
                              R_xlen_t number_types,
                              Configurations&&... configurations) const {
    // TODO: Ideally I'd like to use if constexpr
    // TODO: If I can use if constexpr, I can remove get_nonzero_value below.
    if(has_nonzero_value_v<Varphi>) { // Use faster algorithm in this case
      std::vector<unsigned long long int> count_positive_types(number_types);
      for_each_container([&count_positive_types, this, &point](const auto& current_point) {
        if(count_positive_types[get_type(current_point)] < saturation_) {
          if(Varphi::apply(current_point, point) > 0) {
            count_positive_types[get_type(current_point)] += 1;
          }
        }
      }, std::forward<Configurations>(configurations)...);

      Rcpp::NumericVector dispersion(Rcpp::no_init(number_types));
      for(R_xlen_t i(0); i < number_types; ++i) {
        dispersion[i] = get_nonzero_value<Varphi>() * static_cast<double>(count_positive_types[i]);
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
      Rcpp::NumericVector dispersion(Rcpp::no_init(number_types));
      const auto point_type(get_type(point));
      for(R_xlen_t i(0); i < number_types; ++i) {
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

template<typename Dispersion, typename Lambda, typename Alpha, typename Beta>
class Exponential_family_model: public Dispersion {
public:
  template<typename D, std::enable_if_t<std::is_same<Dispersion, std::remove_reference_t<D>>::value>* = nullptr>
  Exponential_family_model(const Lambda& lambda,
                           const Alpha& alpha,
                           const Beta& beta,
                           Rcpp::List covariates,
                           D&& dispersion):
  Dispersion(std::forward<D>(dispersion)),
    lambda_(lambda),
    alpha_(alpha),
    beta_(beta),
    covariates_(covariates) {}

  template<typename... Args>
  Exponential_family_model(const Lambda& lambda,
                           const Alpha& alpha,
                           const Beta& beta,
                           Rcpp::List covariates,
                           Args&&... args):
    Exponential_family_model(lambda, alpha, beta, covariates, Dispersion(std::forward<Args>(args)...)) {}

  template<typename... Configurations, typename Point>
  double compute_papangelou(const Point& point,
                            R_xlen_t number_types,
                            Configurations&&... configurations) const {
    const auto dispersion(Dispersion::compute(point, number_types, std::forward<Configurations>(configurations)...));

    double inner_product(0);
    const auto point_type(get_type(point));
    for(R_xlen_t i(0); i < number_types; ++i) {
      inner_product += alpha_(point_type, i) * dispersion[i];
    }
    for(decltype(covariates_.size()) i(0); i < covariates_.size(); ++i) {
      inner_product += beta_(point_type, i) * covariates_[i](point);
    }
    return lambda_[point_type] * std::exp(inner_product);
  }

  template<typename Window>
  auto get_dominating_intensity(const Window& window, R_xlen_t number_types) const {
    std::vector<double> dominating_intensity(number_types);
    for(R_xlen_t i(0); i < number_types; ++i) {
      double inner_product(0);
      for(R_xlen_t j(0); j < number_types; ++j) {
        const auto alpha_ij(alpha_(j, i));
        if(alpha_ij > 0) {
          inner_product += alpha_ij * Dispersion::get_maximum(window);
        }
      }
      if(covariates_.size() > 0) {
        inner_product += covariates_.get_maximum_of_dot(beta_(i, Rcpp::_));
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
  Beta beta_;
  Im_list_wrapper covariates_;
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
    return f(Saturated_varphi_model_papangelou<varphi::Exponential>(saturation, radius));
  } else if(model_string == models[1]) {
    return f(Saturated_varphi_model_papangelou<varphi::Square_exponential>(saturation, radius));
  } else if(model_string == models[2]) {
    return f(Saturated_varphi_model_papangelou<varphi::Bump>(saturation, radius));
  } else if(model_string == models[3]) {
    return f(Saturated_varphi_model_papangelou<varphi::Square_bump>(saturation, radius));
  } else if(model_string == models[4]) {
    return f(Saturated_varphi_model_papangelou<varphi::Strauss>(saturation, radius));
  } else {
    Rcpp::stop("Incorrect model entered. A call to show_models() will show you the available choices.\n");
  }
}

template<typename F, typename Lambda>
inline auto call_on_model(Rcpp::CharacterVector model,
                          Rcpp::NumericMatrix alpha,
                          const Lambda& lambda,
                          Rcpp::NumericMatrix beta,
                          Rcpp::List covariates,
                          Rcpp::NumericMatrix radius,
                          unsigned long long int saturation,
                          const F& f) {
  return call_on_dispersion_model(model, radius, saturation, [&alpha, &lambda, &f, &beta, covariates](auto&& varphi) {
    using Model_type = Exponential_family_model<std::remove_cv_t<std::remove_reference_t<decltype(varphi)>>,
                                       Lambda,
                                       Rcpp::NumericMatrix,
                                       Rcpp::NumericMatrix>;
    const Model_type model(lambda, alpha, beta, covariates, std::forward<decltype(varphi)>(varphi));
    return f(std::move(model));
  });
}


} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SATURATED_VARPHI
