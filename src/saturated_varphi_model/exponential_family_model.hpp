#ifndef INCLUDE_PPJSDM_EXPONENTIAL_FAMILY
#define INCLUDE_PPJSDM_EXPONENTIAL_FAMILY

#include <Rcpp.h>

#include "saturated_varphi_model.hpp"
#include "../point/point_manipulation.hpp"
#include "../utility/im_wrapper.hpp"

#include <cmath> // std::exp, std::log, std::fabs, std::isinf
#include <type_traits> // std::remove_cv, std::remove_reference, std::is_same, std::enable_if
#include <utility> // std::forward, std::move
#include <vector> // std::vector

namespace ppjsdm {
namespace detail {

template<typename Point, typename Beta, typename Covariates>
inline auto compute_beta_dot_covariates(const Point& point, const Beta& beta, const Covariates& covariates) {
  double sum(0);
  using size_t = decltype(covariates.size());
  for(size_t i(0); i < covariates.size(); ++i) {
    sum += beta(get_type(point), i) * covariates[i](point);
  }
  return sum;
}

// TODO: Shouldn't this depend on window, i.e. we only compute over covariates in window?
template<typename Beta, typename Covariates>
inline auto compute_beta_dot_covariates_maximum(int type, const Beta& beta, const Covariates& covariates) {
  return covariates.get_maximum_of_dot(beta(type, Rcpp::_));
}

template<typename Beta, typename Covariates>
inline auto compute_beta_dot_covariates_maximum(const Beta& beta, const Covariates& covariates) {
  const auto number_types(beta.nrow());
  std::vector<double> result(number_types);
  if(covariates.size() > 0) {
    for(R_xlen_t i(0); i < number_types; ++i) {
      result[i] = compute_beta_dot_covariates_maximum(i, beta, covariates);
    }
  }
  return result;
}

template<typename Point, typename Alpha, typename Dispersion, typename... Configurations>
inline auto compute_alpha_dot_dispersion(const Point& point,
                                         const Alpha& alpha,
                                         const Dispersion& dispersion,
                                         Configurations&&... configurations) {
  const auto d(dispersion.compute(point, alpha.nrow(), std::forward<Configurations>(configurations)...));
  double sum(0);
  for(R_xlen_t i(0); i < alpha.ncol(); ++i) {
    sum += alpha(get_type(point), i) * d[i];
  }
  return sum;
}

template<typename Point, typename Alpha, typename Dispersion, typename... Configurations>
inline auto compute_positive_alpha_dot_dispersion(const Point& point,
                                                      const Alpha& alpha,
                                                      const Dispersion& dispersion,
                                                      Configurations&&... configurations) {
  const auto d(dispersion.compute(point, alpha.nrow(), std::forward<Configurations>(configurations)...));
  double sum(0);
  for(R_xlen_t i(0); i < alpha.ncol(); ++i) {
    const auto a(alpha(get_type(point), i));
    if(a > 0) {
      sum += a * d[i];
    }
  }
  return sum;
}

template<typename Point, typename Alpha, typename Dispersion, typename... Configurations>
inline auto compute_negative_alpha_dot_dispersion(const Point& point,
                                                  const Alpha& alpha,
                                                  const Dispersion& dispersion,
                                                  Configurations&&... configurations) {
  const auto d(dispersion.compute(point, alpha.nrow(), std::forward<Configurations>(configurations)...));
  double sum(0);
  for(R_xlen_t i(0); i < alpha.ncol(); ++i) {
    const auto a(alpha(get_type(point), i));
    if(a < 0) {
      sum += a * d[i];
    }
  }
  return sum;
}

template<typename Point, typename Alpha>
inline bool is_alpha_non_negative(const Point& point,
                                  const Alpha& alpha) {
  for(R_xlen_t i(0); i < alpha.ncol(); ++i) {
    if(alpha(get_type(point), i) < 0) {
      return false;
    }
  }
  return true;
}

template<typename Point, typename Alpha>
inline bool is_alpha_non_positive(const Point& point,
                                  const Alpha& alpha) {
  for(R_xlen_t i(0); i < alpha.ncol(); ++i) {
    if(alpha(get_type(point), i) > 0) {
      return false;
    }
  }
  return true;
}

template<typename Alpha>
inline auto compute_alpha_dot_dispersion_maximum(const Alpha& alpha,
                                                 double dispersion_maximum,
                                                 int type) {
  double sum(0);
  for(R_xlen_t j(0); j < alpha.ncol(); ++j) {
    const auto a(alpha(type, j));
    if(a > 0) {
      sum += a;
    }
  }
  return sum * dispersion_maximum;
}

template<typename Alpha>
inline auto compute_alpha_dot_dispersion_maximum(const Alpha& alpha, double dispersion_maximum) {
  const auto number_types(alpha.nrow());
  std::vector<double> result(number_types);
  for(R_xlen_t i(0); i < number_types; ++i) {
    result[i] = compute_alpha_dot_dispersion_maximum(alpha, dispersion_maximum, i);
  }
  return result;
}

} // namespace detail

template<typename Dispersion, typename Lambda>
class Exponential_family_model {
public:
  template<typename D, std::enable_if_t<std::is_same<Dispersion, std::remove_reference_t<D>>::value>* = nullptr>
  Exponential_family_model(const Lambda& lambda,
                           D&& dispersion,
                           Rcpp::NumericMatrix alpha,
                           Rcpp::NumericMatrix beta,
                           Rcpp::List covariates):
    dispersion_(std::forward<D>(dispersion)),
    lambda_(lambda),
    alpha_(alpha),
    beta_(beta),
    covariates_(covariates) {}

  template<typename Point, typename... Configurations>
  double compute_papangelou(const Point& point, Configurations&&... configurations) const {
    double alpha_dispersion(detail::compute_alpha_dot_dispersion(point, alpha_, dispersion_, std::forward<Configurations>(configurations)...));
    double beta_covariates(detail::compute_beta_dot_covariates(point, beta_, covariates_));
    return lambda_[get_type(point)] * std::exp(alpha_dispersion + beta_covariates);
  }

protected:
  Dispersion dispersion_;
  Lambda lambda_;
  Rcpp::NumericMatrix alpha_;
  Rcpp::NumericMatrix beta_;
  Im_list_wrapper covariates_;
};

template<typename Window, typename Dispersion, typename Lambda>
class Exponential_family_model_over_window: public Exponential_family_model<Dispersion, Lambda> {
private:
  using Model = Exponential_family_model<Dispersion, Lambda>;
public:
  template<typename D>
  Exponential_family_model_over_window(const Window& window,
                                       const Lambda& lambda,
                                       D&& dispersion,
                                       Rcpp::NumericMatrix alpha,
                                       Rcpp::NumericMatrix beta,
                                       Rcpp::List covariates):
    Model(lambda, std::forward<D>(dispersion), alpha, beta, covariates),
    window_(window),
    beta_dot_covariates_maximum_(detail::compute_beta_dot_covariates_maximum(beta, Model::covariates_)),
    alpha_dot_dispersion_maximum_(detail::compute_alpha_dot_dispersion_maximum(alpha, dispersion.get_maximum(window))) {}

  template<typename Point>
  auto get_log_normalized_bounding_intensity(const Point& point) const {
    double beta_covariates_maximum(beta_dot_covariates_maximum_[get_type(point)]);
    double beta_covariates(detail::compute_beta_dot_covariates(point, Model::beta_, Model::covariates_));
    return beta_covariates - beta_covariates_maximum;
  }

  // TODO: This should somehow be restricted to window.
  auto get_integral() const {
    double integral(0);
    const auto number_types(Model::beta_.nrow());
    for(R_xlen_t i(0); i < number_types; ++i) {
      integral += Model::covariates_.get_integral_of_dot([](double x) { return std::exp(x); }, Model::beta_(i, Rcpp::_));
    }
    return integral / static_cast<double>(number_types);
  }

  auto sample_point_from_bounding_intensity(R_xlen_t type) const {
    while(true) {
      const auto sample(window_.sample(type));
      if(exp_rand() + get_log_normalized_bounding_intensity(sample) >= 0) {
        return sample;
      }
    }
  }

  auto get_upper_bound() const {
    const auto number_types(Model::lambda_.size());
    std::vector<double> upper_bound(number_types);
    using size_t = decltype(Model::lambda_.size());
    for(size_t i(0); i < number_types; ++i) {
      const auto value(Model::lambda_[i] * std::exp(alpha_dot_dispersion_maximum_[i] + beta_dot_covariates_maximum_[i]));
      if(std::isinf(value)) {
        Rcpp::stop("Infinite value obtained as the bound to the Papangelou intensity.");
      }
      upper_bound[i] = value;
    }
    return upper_bound;
  }

  template<typename Point, typename... Configurations>
  double compute_normalised_papangelou(const Point& point, Configurations&&... configurations) const {
    double alpha_dispersion_maximum(detail::compute_alpha_dot_dispersion_maximum(Model::alpha_, Model::dispersion_.get_maximum(window_), get_type(point)));
    double alpha_dispersion(detail::compute_alpha_dot_dispersion(point,Model:: alpha_, Model::dispersion_, std::forward<Configurations>(configurations)...));
    return std::exp(alpha_dispersion - alpha_dispersion_maximum);
  }

  template<typename Point, typename L, typename LComplement>
  double compute_log_alpha_max(const Point& point, const L& l, const LComplement& l_complement) const {
    double alpha_dispersion_maximum(detail::compute_alpha_dot_dispersion_maximum(Model::alpha_, Model::dispersion_.get_maximum(window_), get_type(point)));
    double alpha_dispersion;
    if(detail::is_alpha_non_negative(point, Model::alpha_)) {
      alpha_dispersion = detail::compute_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l, l_complement);
    } else if(detail::is_alpha_non_positive(point, Model::alpha_)) {
      alpha_dispersion = detail::compute_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l);
    } else {
      const auto positive_alpha(detail::compute_positive_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l, l_complement));
      const auto negative_alpha(detail::compute_negative_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l));
      alpha_dispersion = positive_alpha + negative_alpha;
    }
    return alpha_dispersion - alpha_dispersion_maximum;
  }

  // TODO: Factorise this and above.
  template<typename Point, typename L, typename LComplement>
  double compute_log_alpha_min(const Point& point, const L& l, const LComplement& l_complement) const {
    double alpha_dispersion_maximum(detail::compute_alpha_dot_dispersion_maximum(Model::alpha_, Model::dispersion_.get_maximum(window_), get_type(point)));
    double alpha_dispersion;
    if(detail::is_alpha_non_negative(point, Model::alpha_)) {
      alpha_dispersion = detail::compute_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l);
    } else if(detail::is_alpha_non_positive(point, Model::alpha_)) {
      alpha_dispersion = detail::compute_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l, l_complement);
    } else {
      const auto positive_alpha(detail::compute_positive_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l));
      const auto negative_alpha(detail::compute_negative_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l, l_complement));
      alpha_dispersion = positive_alpha + negative_alpha;
    }
    return alpha_dispersion - alpha_dispersion_maximum;
  }

  double compute_log_alpha_min_lower_bound(R_xlen_t type) const {
    double sum(0);
    for(R_xlen_t j(0); j < Model::alpha_.ncol(); ++j) {
      const auto a(Model::alpha_(type, j));
      if(a > 0) {
        sum -= a;
      } else {
        sum += a;
      }
    }
    return sum * Model::dispersion_.get_maximum(window_);
  }

  const auto& get_window() const {
    return window_;
  }

private:
  Window window_;
  std::vector<double> beta_dot_covariates_maximum_;
  std::vector<double> alpha_dot_dispersion_maximum_;
};

template<typename F, typename Lambda, typename... Args>
inline auto call_on_model(Rcpp::CharacterVector model,
                          const Lambda& lambda,
                          Rcpp::NumericMatrix radius,
                          unsigned long long int saturation,
                          const F& f,
                          Args... args) {
  return call_on_dispersion_model(model, radius, saturation, [&lambda, &f, args...](auto&& varphi) {
    using Model_type = Exponential_family_model<std::remove_cv_t<std::remove_reference_t<decltype(varphi)>>, Lambda>;
    return f(Model_type(lambda, std::forward<decltype(varphi)>(varphi), args...));
  });
}

template<typename Window, typename F, typename Lambda, typename... Args>
inline auto call_on_model(const Window& window,
                          Rcpp::CharacterVector model,
                          const Lambda& lambda,
                          Rcpp::NumericMatrix radius,
                          unsigned long long int saturation,
                          const F& f,
                          Args... args) {
  return call_on_dispersion_model(model, radius, saturation, [&window, &lambda, &f, args...](auto&& varphi) {
    using Model_type = Exponential_family_model_over_window<Window,
                                                            std::remove_cv_t<std::remove_reference_t<decltype(varphi)>>,
                                                            Lambda>;
    return f(Model_type(window, lambda, std::forward<decltype(varphi)>(varphi), args...));
  });
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_EXPONENTIAL_FAMILY
