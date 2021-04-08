#ifndef INCLUDE_PPJSDM_EXPONENTIAL_FAMILY
#define INCLUDE_PPJSDM_EXPONENTIAL_FAMILY

#include <Rcpp.h>

#include "saturated_model.hpp"
#include "../configuration/configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"
#include "../simulation/inhomogeneous_ppp.hpp"
#include "../utility/im_wrapper.hpp"
#include "../utility/approximate_draw.hpp"

#include <algorithm> // std::transform
#include <cmath> // std::exp, std::isinf
#include <functional> // std::plus
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

template<typename Point, typename Alpha>
inline bool is_alpha_zero(const Point& point,
                          const Alpha& alpha) {
  for(R_xlen_t i(0); i < alpha.ncol(); ++i) {
    if(alpha(get_type(point), i) != 0) {
      return false;
    }
  }
  return true;
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

template<typename Point, typename Alpha, typename Dispersion, typename... Configurations>
inline auto compute_alpha_dot_dispersion(const Point& point,
                                         const Alpha& alpha,
                                         const Dispersion& dispersion,
                                         Configurations&... configurations) {
  if(is_alpha_zero(point, alpha)) {
    return 0.;
  }
  const auto dispersion_value(compute_dispersion(dispersion, point, alpha.nrow(), configurations...));
  double sum(0);
  for(R_xlen_t i(0); i < alpha.ncol(); ++i) {
    if(dispersion_value[i] != 0) { // This ensures that \infty \times 0 = 0.
      sum += alpha(get_type(point), i) * dispersion_value[i];
    }
  }
  return sum;
}

template<typename Point, typename Alpha, typename Dispersion, typename... Configurations>
inline auto compute_positive_alpha_dot_dispersion(const Point& point,
                                                  const Alpha& alpha,
                                                  const Dispersion& dispersion,
                                                  Configurations&... configurations) {
  if(is_alpha_zero(point, alpha)) {
    return 0.;
  }
  const auto dispersion_value(compute_dispersion(dispersion, point, alpha.nrow(), configurations...));
  double sum(0);
  for(R_xlen_t i(0); i < alpha.ncol(); ++i) {
    const auto a(alpha(get_type(point), i));
    if(a > 0) {
      sum += a * dispersion_value[i];
    }
  }
  return sum;
}

template<typename Point, typename Alpha, typename Dispersion, typename... Configurations>
inline auto compute_negative_alpha_dot_dispersion(const Point& point,
                                                  const Alpha& alpha,
                                                  const Dispersion& dispersion,
                                                  Configurations&... configurations) {
  if(is_alpha_zero(point, alpha)) {
    return 0.;
  }
  const auto dispersion_value(compute_dispersion(dispersion, point, alpha.nrow(), configurations...));
  double sum(0);
  for(R_xlen_t i(0); i < alpha.ncol(); ++i) {
    const auto a(alpha(get_type(point), i));
    if(a < 0) {
      sum += a * dispersion_value[i];
    }
  }
  return sum;
}

template<typename Alpha>
inline auto compute_alpha_dot_dispersion_maximum(const Alpha& alpha,
                                                 double dispersion_maximum,
                                                 int type) {
  double sum(0.);
  for(R_xlen_t j(0); j < alpha.ncol(); ++j) {
    const auto a(alpha(type, j));
    if(a > 0.) {
      sum += a;
    }
  }
  // The condition below ensures that 0 * Inf = 0.
  if(sum > 0.) {
    return sum * dispersion_maximum;
  } else {
    return 0.;
  }
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

template<typename Lambda>
class Truncated_exponential_family_model {
public:
  Truncated_exponential_family_model(const Lambda& beta0,
                                     Rcpp::CharacterVector model,
                                     Rcpp::CharacterVector medium_range_model,
                                     Rcpp::NumericMatrix alpha,
                                     Rcpp::NumericMatrix beta,
                                     Rcpp::NumericMatrix gamma,
                                     Rcpp::List covariates,
                                     Rcpp::NumericMatrix short_range,
                                     Rcpp::NumericMatrix medium_range,
                                     Rcpp::NumericMatrix long_range,
                                     unsigned long long int saturation):
    dispersion_(Saturated_model(model, short_range, saturation)),
    medium_range_dispersion_(Saturated_model(medium_range_model, medium_range, long_range, saturation)),
    beta0_(beta0),
    alpha_(alpha),
    beta_(beta),
    gamma_(gamma),
    covariates_(covariates) {}

  template<typename Point, typename Configuration>
  double compute_log_papangelou(const Point& point, const Configuration& configuration) const {
    double alpha_dispersion(detail::compute_alpha_dot_dispersion(point, alpha_, dispersion_, configuration));
    double gamma_dispersion(detail::compute_alpha_dot_dispersion(point, gamma_, medium_range_dispersion_, configuration));
    double beta_covariates(detail::compute_beta_dot_covariates(point, beta_, covariates_));
    return beta0_[get_type(point)] + beta_covariates + alpha_dispersion + gamma_dispersion;
  }

  template<typename Point, typename Configuration>
  double compute_papangelou(const Point& point, const Configuration& configuration) const {
    return std::exp(compute_log_papangelou(point, configuration));
  }

  auto get_number_types() const {
    return beta0_.size();
  }

protected:
  Saturated_model dispersion_;
  Saturated_model medium_range_dispersion_;
  Lambda beta0_;
  Rcpp::NumericMatrix alpha_;
  Rcpp::NumericMatrix beta_;
  Rcpp::NumericMatrix gamma_;
  Im_list_wrapper covariates_;
};

template<typename Lambda>
class Truncated_exponential_family_model_over_window: public Truncated_exponential_family_model<Lambda> {
private:
  using Model = Truncated_exponential_family_model<Lambda>;

  template<typename D>
  auto get_dispersion_maximum(const D& dispersion) const {
    if(dispersion.is_nonincreasing()) {
      return 6 * dispersion.get_maximum();
    } else {
      return std::numeric_limits<double>::infinity();
    }
  }
public:
  Truncated_exponential_family_model_over_window(const Window& window,
                                                 const Lambda& beta0,
                                                 Rcpp::CharacterVector model,
                                                 Rcpp::CharacterVector medium_range_model,
                                                 Rcpp::NumericMatrix alpha,
                                                 Rcpp::NumericMatrix beta,
                                                 Rcpp::NumericMatrix gamma,
                                                 Rcpp::List covariates,
                                                 Rcpp::NumericMatrix short_range,
                                                 Rcpp::NumericMatrix medium_range,
                                                 Rcpp::NumericMatrix long_range,
                                                 unsigned long long int saturation):
    Model(beta0, model, medium_range_model, alpha, beta, gamma, covariates, short_range, medium_range, long_range, saturation),
    window_(window),
    beta_dot_covariates_maximum_(detail::compute_beta_dot_covariates_maximum(beta, Model::covariates_)),
    dot_dispersion_maximum_(detail::compute_alpha_dot_dispersion_maximum(alpha, get_dispersion_maximum(Model::dispersion_))) {
    const auto gamma_dot_dispersion_maximum(detail::compute_alpha_dot_dispersion_maximum(gamma, get_dispersion_maximum(Model::medium_range_dispersion_)));
    std::transform(dot_dispersion_maximum_.begin(), dot_dispersion_maximum_.end(), gamma_dot_dispersion_maximum.begin(),
                   dot_dispersion_maximum_.begin(), std::plus<double>());
  }

  template<typename Point>
  auto get_log_normalized_bounding_intensity(const Point& point) const {
    double beta_covariates_maximum(beta_dot_covariates_maximum_[get_type(point)]);
    double beta_covariates(detail::compute_beta_dot_covariates(point, Model::beta_, Model::covariates_));
    return beta_covariates - beta_covariates_maximum;
  }

  template<typename Point>
  auto get_log_approximate_ppp_intensity(const Point& point) const {
    return get_log_normalized_bounding_intensity(point);
  }

  // TODO: This should somehow be restricted to window.
  auto get_integral() const {
    double integral(0);
    const auto number_types(Model::beta_.nrow());
    for(R_xlen_t i(0); i < number_types; ++i) {
      integral += Model::covariates_.get_integral_of_dot(window_, [i, this](double x) {
        return std::exp(x + Model::beta0_[i] + dot_dispersion_maximum_[i]);
      }, Model::beta_(i, Rcpp::_));
    }
    return integral;
  }

  auto sample_point_from_bounding_intensity() const {
    // Sample type proportionally to the exp(beta0_i).
    const auto random_type(Rcpp::sample(Model::beta0_.size(), 1, false, Rcpp::sugar::probs_t(Rcpp::exp(Model::beta0_)), false)[0]);
    while(true) {
      const auto sample(window_.sample(random_type));
      if(exp_rand() + get_log_normalized_bounding_intensity(sample) >= 0) {
        return sample;
      }
    }
  }

  auto get_upper_bound_approximate_ppp_intensity() const {
    const auto number_types(Model::beta0_.size());
    std::vector<double> upper_bound(number_types);
    using size_t = decltype(Model::beta0_.size());
    for(size_t i(0); i < number_types; ++i) {
      const auto value(std::exp(Model::beta0_[i] + beta_dot_covariates_maximum_[i]));
      if(std::isinf(value)) {
        Rcpp::stop("Infinite value obtained as the bound to the approximate PPP intensity.");
      }
      upper_bound[i] = value;
    }
    return upper_bound;
  }

  auto get_upper_bound() const {
    const auto number_types(Model::beta0_.size());
    std::vector<double> upper_bound(number_types);
    using size_t = decltype(Model::beta0_.size());
    for(size_t i(0); i < number_types; ++i) {
      const auto value(std::exp(Model::beta0_[i] + dot_dispersion_maximum_[i] + beta_dot_covariates_maximum_[i]));
      if(std::isinf(value)) {
        Rcpp::Rcout << Model::beta0_[i] << '\n';
        Rcpp::Rcout << dot_dispersion_maximum_[i] << '\n';
        Rcpp::Rcout << beta_dot_covariates_maximum_[i] << '\n';
        Rcpp::stop("Infinite value obtained as the bound to the Papangelou intensity.");
      }
      upper_bound[i] = value;
    }
    return upper_bound;
  }

  // TODO: Factorise and make this function more readable.
  template<typename Point, typename Configuration>
  void add_to_L_or_U(double exp_mark, const Point& point, Configuration& l, Configuration& l_complement) const {
    // TODO: Write properly what this function is doing.
    double dot_dispersion_maximum(dot_dispersion_maximum_[get_type(point)]);
    if(detail::is_alpha_non_negative(point, Model::alpha_) && detail::is_alpha_non_negative(point, Model::gamma_)) {
      const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l, l_complement));
      const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, Model::medium_range_dispersion_, l, l_complement));
      const auto log_alpha(alpha_dispersion + gamma_dispersion - dot_dispersion_maximum);
      if(log_alpha > 0) {
        Rcpp::stop("Bad upper-bound in the computation of alpha in the CFTP algorithm.");
      }
      if(log_alpha + exp_mark > 0) {
        const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l));
        const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, Model::medium_range_dispersion_, l));
        const auto log_alpha(alpha_dispersion + gamma_dispersion - dot_dispersion_maximum);
        if(log_alpha > 0) {
          Rcpp::stop("Bad upper-bound in the computation of alpha in the CFTP algorithm.");
        }
        add_point(log_alpha + exp_mark > 0 ? l : l_complement, point);
      }
    } else if(detail::is_alpha_non_positive(point, Model::alpha_) && detail::is_alpha_non_positive(point, Model::gamma_)) {
      const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l));
      const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, Model::medium_range_dispersion_, l));
      const auto log_alpha(alpha_dispersion + gamma_dispersion - dot_dispersion_maximum);
      if(log_alpha > 0) {
        Rcpp::stop("Bad upper-bound in the computation of alpha in the CFTP algorithm.");
      }
      if(log_alpha + exp_mark > 0) {
        const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l, l_complement));
        const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, Model::medium_range_dispersion_, l, l_complement));
        const auto log_alpha(alpha_dispersion + gamma_dispersion - dot_dispersion_maximum);
        if(log_alpha > 0) {
          Rcpp::stop("Bad upper-bound in the computation of alpha in the CFTP algorithm.");
        }
        add_point(log_alpha + exp_mark > 0 ? l : l_complement, point);
      }
    } else {
      const auto positive_alpha(detail::compute_positive_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l, l_complement));
      const auto negative_alpha(detail::compute_negative_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l));
      const auto positive_gamma(detail::compute_positive_alpha_dot_dispersion(point, Model::gamma_, Model::medium_range_dispersion_, l, l_complement));
      const auto negative_gamma(detail::compute_negative_alpha_dot_dispersion(point, Model::gamma_, Model::medium_range_dispersion_, l));
      const auto log_alpha(positive_alpha + negative_alpha + positive_gamma + negative_gamma - dot_dispersion_maximum);
      if(log_alpha > 0) {
        Rcpp::stop("Bad upper-bound in the computation of alpha in the CFTP algorithm.");
      }
      if(positive_alpha + negative_alpha + positive_gamma + negative_gamma - dot_dispersion_maximum + exp_mark > 0) {
        const auto positive_alpha(detail::compute_positive_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l));
        const auto negative_alpha(detail::compute_negative_alpha_dot_dispersion(point, Model::alpha_, Model::dispersion_, l, l_complement));
        const auto positive_gamma(detail::compute_positive_alpha_dot_dispersion(point, Model::gamma_, Model::medium_range_dispersion_, l));
        const auto negative_gamma(detail::compute_negative_alpha_dot_dispersion(point, Model::gamma_, Model::medium_range_dispersion_, l, l_complement));
        const auto log_alpha(positive_alpha + negative_alpha + positive_gamma + negative_gamma - dot_dispersion_maximum);
        if(log_alpha > 0) {
          Rcpp::stop("Bad upper-bound in the computation of alpha in the CFTP algorithm.");
        }
        add_point(log_alpha + exp_mark > 0 ? l : l_complement, point);
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

    return sum_alpha * get_dispersion_maximum(Model::dispersion_) + sum_gamma * get_dispersion_maximum(Model::medium_range_dispersion_);
  }

  const auto& get_window() const {
    return window_;
  }

private:
  Window window_;
  std::vector<double> beta_dot_covariates_maximum_;
  std::vector<double> dot_dispersion_maximum_;
};

namespace detail {

template<typename Configuration, typename Lambda>
struct approximate_draw_helper<Configuration, Truncated_exponential_family_model_over_window<Lambda>> {
  static auto get(const Truncated_exponential_family_model_over_window<Lambda>& model) {
    const auto configuration(simulate_inhomogeneous_ppp<Configuration>(model.get_window(),
                                                                       [&model](const auto& point) {
                                                                         return model.get_log_approximate_ppp_intensity(point);
                                                                       },
                                                                       model.get_upper_bound_approximate_ppp_intensity(),
                                                                       model.get_number_types()));
    return configuration;
  }
};

} // namespace detail

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_EXPONENTIAL_FAMILY
