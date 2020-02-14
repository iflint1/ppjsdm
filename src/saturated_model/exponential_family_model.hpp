#ifndef INCLUDE_PPJSDM_EXPONENTIAL_FAMILY
#define INCLUDE_PPJSDM_EXPONENTIAL_FAMILY

#include <Rcpp.h>

#include "saturated_model.hpp"
#include "../configuration/configuration_manipulation.hpp"
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

template<typename Point, typename Alpha, typename Dispersion>
inline auto compute_alpha_dot_dispersion(const Point& point,
                                         const Alpha& alpha,
                                         const Dispersion& dispersion) {
  double sum(0);
  for(R_xlen_t i(0); i < alpha.ncol(); ++i) {
    sum += alpha(get_type(point), i) * dispersion[i];
  }
  return sum;
}

template<typename Point, typename Alpha, typename Dispersion>
inline auto compute_positive_alpha_dot_dispersion(const Point& point,
                                                  const Alpha& alpha,
                                                  const Dispersion& dispersion) {
  double sum(0);
  for(R_xlen_t i(0); i < alpha.ncol(); ++i) {
    const auto a(alpha(get_type(point), i));
    if(a > 0) {
      sum += a * dispersion[i];
    }
  }
  return sum;
}

template<typename Point, typename Alpha, typename Dispersion>
inline auto compute_negative_alpha_dot_dispersion(const Point& point,
                                                  const Alpha& alpha,
                                                  const Dispersion& dispersion) {
  double sum(0);
  for(R_xlen_t i(0); i < alpha.ncol(); ++i) {
    const auto a(alpha(get_type(point), i));
    if(a < 0) {
      sum += a * dispersion[i];
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

template<typename Dispersion, typename MediumRangeDispersion, typename Lambda>
class Exponential_family_model {
public:
  template<typename D, typename E, std::enable_if_t<std::is_same<Dispersion, std::remove_reference_t<D>>::value && std::is_same<MediumRangeDispersion, std::remove_reference_t<E>>::value>* = nullptr>
  Exponential_family_model(const Lambda& lambda,
                           D&& dispersion,
                           E&& medium_range_dispersion,
                           Rcpp::NumericMatrix alpha,
                           Rcpp::NumericMatrix beta,
                           Rcpp::NumericMatrix gamma,
                           Rcpp::List covariates):
    dispersion_(std::forward<D>(dispersion)),
    medium_range_dispersion_(std::forward<E>(medium_range_dispersion)),
    lambda_(lambda),
    alpha_(alpha),
    beta_(beta),
    gamma_(gamma),
    covariates_(covariates) {}

  template<typename Point, typename Configuration>
  double compute_papangelou(const Point& point, const Configuration& configuration) const {
    const auto dalpha(dispersion_.compute(point, alpha_.nrow(), configuration));
    double alpha_dispersion(detail::compute_alpha_dot_dispersion(point, alpha_, dalpha));
    const auto dgamma(medium_range_dispersion_.compute(point, alpha_.nrow(), configuration));
    double gamma_dispersion(detail::compute_alpha_dot_dispersion(point, gamma_, dgamma));
    double beta_covariates(detail::compute_beta_dot_covariates(point, beta_, covariates_));
    return lambda_[get_type(point)] * std::exp(alpha_dispersion + beta_covariates + gamma_dispersion);
  }

protected:
  Dispersion dispersion_;
  MediumRangeDispersion medium_range_dispersion_;
  Lambda lambda_;
  Rcpp::NumericMatrix alpha_;
  Rcpp::NumericMatrix beta_;
  Rcpp::NumericMatrix gamma_;
  Im_list_wrapper covariates_;
};

template<typename Window, typename Dispersion, typename MediumRangeDispersion, typename Lambda>
class Exponential_family_model_over_window: public Exponential_family_model<Dispersion, MediumRangeDispersion, Lambda> {
private:
  using Model = Exponential_family_model<Dispersion, MediumRangeDispersion, Lambda>;
public:
  template<typename D, typename E>
  Exponential_family_model_over_window(const Window& window,
                                       const Lambda& lambda,
                                       D&& dispersion,
                                       E&& medium_range_dispersion,
                                       Rcpp::NumericMatrix alpha,
                                       Rcpp::NumericMatrix beta,
                                       Rcpp::NumericMatrix gamma,
                                       Rcpp::List covariates):
    Model(lambda, std::forward<D>(dispersion), std::forward<E>(medium_range_dispersion), alpha, beta, gamma, covariates),
    window_(window),
    beta_dot_covariates_maximum_(detail::compute_beta_dot_covariates_maximum(beta, Model::covariates_)),
    alpha_dot_dispersion_maximum_(detail::compute_alpha_dot_dispersion_maximum(alpha, dispersion.get_maximum(window))),
    gamma_dot_dispersion_maximum_(detail::compute_alpha_dot_dispersion_maximum(gamma, medium_range_dispersion.get_maximum(window))){}

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
      const auto value(Model::lambda_[i] * std::exp(alpha_dot_dispersion_maximum_[i] + gamma_dot_dispersion_maximum_[i] + beta_dot_covariates_maximum_[i]));
      if(std::isinf(value)) {
        Rcpp::stop("Infinite value obtained as the bound to the Papangelou intensity.");
      }
      upper_bound[i] = value;
    }
    return upper_bound;
  }

  // TODO: Factorise and make this function more readable.
  // TODO: Does not take into account gamma dispersion
  template<typename Point, typename Configuration>
  void add_to_L_or_U(double exp_mark, const Point& point, Configuration& l, Configuration& l_complement) const {
    double alpha_dispersion_maximum(detail::compute_alpha_dot_dispersion_maximum(Model::alpha_, Model::dispersion_.get_maximum(window_), get_type(point)));
    double gamma_dispersion_maximum(detail::compute_alpha_dot_dispersion_maximum(Model::gamma_, Model::medium_range_dispersion_.get_maximum(window_), get_type(point)));
    auto compute_alpha_dispersion(Model::dispersion_.get_compute_dispersion_object(point, Model::alpha_.nrow()));
    auto compute_gamma_dispersion(Model::medium_range_dispersion_.get_compute_dispersion_object(point, Model::gamma_.nrow()));
    compute_alpha_dispersion.add_configuration(l);
    compute_gamma_dispersion.add_configuration(l);
    if(detail::is_alpha_non_negative(point, Model::alpha_) && detail::is_alpha_non_negative(point, Model::gamma_)) {
      const auto save_alpha(compute_alpha_dispersion.get_state());
      const auto save_gamma(compute_gamma_dispersion.get_state());
      compute_alpha_dispersion.add_configuration(l_complement);
      compute_gamma_dispersion.add_configuration(l_complement);
      const auto alpha_dispersion_u(compute_alpha_dispersion.compute());
      const auto gamma_dispersion_u(compute_gamma_dispersion.compute());
      const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_u));
      const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_u));
      if(alpha_dispersion + gamma_dispersion - alpha_dispersion_maximum - gamma_dispersion_maximum + exp_mark > 0) {
        const auto alpha_dispersion_l(compute_alpha_dispersion.compute_from_state(save_alpha));
        const auto gamma_dispersion_l(compute_gamma_dispersion.compute_from_state(save_gamma));
        const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_l));
        const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_l));
        add_point(alpha_dispersion + gamma_dispersion - alpha_dispersion_maximum - gamma_dispersion_maximum + exp_mark > 0 ? l : l_complement, point);
      }
    } else if(detail::is_alpha_non_positive(point, Model::alpha_) && detail::is_alpha_non_positive(point, Model::gamma_)) {
      const auto alpha_dispersion_l(compute_alpha_dispersion.compute());
      const auto gamma_dispersion_l(compute_gamma_dispersion.compute());
      const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_l));
      const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_l));
      if(alpha_dispersion + gamma_dispersion - alpha_dispersion_maximum - gamma_dispersion_maximum + exp_mark > 0) {
        compute_alpha_dispersion.add_configuration(l_complement);
        compute_gamma_dispersion.add_configuration(l_complement);
        const auto alpha_dispersion_u(compute_alpha_dispersion.compute());
        const auto gamma_dispersion_u(compute_gamma_dispersion.compute());
        const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_u));
        const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_u));
        add_point(alpha_dispersion + gamma_dispersion - alpha_dispersion_maximum - gamma_dispersion_maximum + exp_mark > 0 ? l : l_complement, point);
      }
    } else {
      const auto alpha_dispersion_l(compute_alpha_dispersion.compute());
      const auto gamma_dispersion_l(compute_gamma_dispersion.compute());
      compute_alpha_dispersion.add_configuration(l_complement);
      compute_gamma_dispersion.add_configuration(l_complement);
      const auto alpha_dispersion_u(compute_alpha_dispersion.compute());
      const auto gamma_dispersion_u(compute_gamma_dispersion.compute());
      const auto positive_alpha(detail::compute_positive_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_u));
      const auto negative_alpha(detail::compute_negative_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_l));
      const auto positive_gamma(detail::compute_positive_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_u));
      const auto negative_gamma(detail::compute_negative_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_l));
      if(positive_alpha + negative_alpha + positive_gamma + negative_gamma - alpha_dispersion_maximum - gamma_dispersion_maximum + exp_mark > 0) {
        const auto positive_alpha(detail::compute_positive_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_l));
        const auto negative_alpha(detail::compute_negative_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_u));
        const auto positive_gamma(detail::compute_positive_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_l));
        const auto negative_gamma(detail::compute_negative_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_u));
        add_point(positive_alpha + negative_alpha + positive_gamma + negative_gamma - alpha_dispersion_maximum - gamma_dispersion_maximum + exp_mark > 0 ? l : l_complement, point);
      }
    }
  }

  double compute_log_alpha_min_lower_bound(R_xlen_t type) const {
    double sum_alpha(0);
    for(R_xlen_t j(0); j < Model::alpha_.ncol(); ++j) {
      const auto a(Model::alpha_(type, j));
      if(a > 0) {
        sum_alpha -= a;
      } else {
        sum_alpha += a;
      }
    }
    double sum_gamma(0);
    for(R_xlen_t j(0); j < Model::gamma_.ncol(); ++j) {
      const auto g(Model::gamma_(type, j));
      if(g > 0) {
        sum_gamma -= g;
      } else {
        sum_gamma += g;
      }
    }

    return sum_alpha * Model::dispersion_.get_maximum(window_)
      + sum_gamma * Model::medium_range_dispersion_.get_maximum(window_);
  }

  const auto& get_window() const {
    return window_;
  }

private:
  Window window_;
  std::vector<double> beta_dot_covariates_maximum_;
  std::vector<double> alpha_dot_dispersion_maximum_;
  std::vector<double> gamma_dot_dispersion_maximum_;
};

template<typename F, typename Lambda, typename... Args>
inline auto call_on_model(Rcpp::CharacterVector model,
                          Rcpp::CharacterVector medium_range_model,
                          const Lambda& lambda,
                          Rcpp::NumericMatrix short_range,
                          Rcpp::NumericMatrix medium_range,
                          Rcpp::NumericMatrix long_range,
                          unsigned long long int saturation,
                          const F& f,
                          Args... args) {
  return call_on_dispersion_model(model, short_range, saturation, [medium_range_model, medium_range, long_range, saturation, &lambda, &f, args...](auto&& varphi) {
    return call_on_medium_range_dispersion_model(medium_range_model, medium_range, long_range, saturation, [&varphi, &lambda, &f, args...](auto&& psi) {

    using Model_type = Exponential_family_model<std::remove_cv_t<std::remove_reference_t<decltype(varphi)>>,
                                                std::remove_cv_t<std::remove_reference_t<decltype(psi)>>,
                                                Lambda>;
    return f(Model_type(lambda, std::forward<decltype(varphi)>(varphi), std::forward<decltype(psi)>(psi), args...));
    });
  });
}

template<typename Window, typename F, typename Lambda, typename... Args>
inline auto call_on_model(const Window& window,
                          Rcpp::CharacterVector model,
                          Rcpp::CharacterVector medium_range_model,
                          const Lambda& lambda,
                          Rcpp::NumericMatrix short_range,
                          Rcpp::NumericMatrix medium_range,
                          Rcpp::NumericMatrix long_range,
                          unsigned long long int saturation,
                          const F& f,
                          Args... args) {
  return call_on_dispersion_model(model, short_range, saturation, [medium_range_model, medium_range, long_range, saturation, &window, &lambda, &f, args...](auto&& varphi) {
    return call_on_medium_range_dispersion_model(medium_range_model, medium_range, long_range, saturation, [&varphi, &window, &lambda, &f, args...](auto&& psi) {
      using Model_type = Exponential_family_model_over_window<Window,
                                                              std::remove_cv_t<std::remove_reference_t<decltype(varphi)>>,
                                                              std::remove_cv_t<std::remove_reference_t<decltype(psi)>>,
                                                              Lambda>;
      return f(Model_type(window, lambda, std::forward<decltype(varphi)>(varphi), std::forward<decltype(psi)>(psi), args...));
    });
  });
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_EXPONENTIAL_FAMILY