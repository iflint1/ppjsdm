#ifndef INCLUDE_PPJSDM_EXPONENTIAL_FAMILY
#define INCLUDE_PPJSDM_EXPONENTIAL_FAMILY

#include <Rcpp.h>

#include "saturated_model.hpp"
#include "../configuration/configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"
#include "../utility/im_wrapper.hpp"

#include <algorithm> // std::transform
#include <cmath> // std::exp, std::log, std::fabs, std::isinf
#include <functional> // std::plus
#include <memory> // std::shared_ptr
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
    if(dispersion[i] != 0) { // This ensures that \infty \times 0 = 0.
      sum += alpha(get_type(point), i) * dispersion[i];
    }
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

// template<typename Dispersion, typename MediumRangeDispersion, typename Lambda>
// class Exponential_family_model {
// public:
//   template<typename D, typename E, std::enable_if_t<std::is_same<Dispersion, std::remove_reference_t<D>>::value && std::is_same<MediumRangeDispersion, std::remove_reference_t<E>>::value>* = nullptr>
//   Exponential_family_model(const Lambda& lambda,
//                            D&& dispersion,
//                            E&& medium_range_dispersion,
//                            Rcpp::NumericMatrix alpha,
//                            Rcpp::NumericMatrix beta,
//                            Rcpp::NumericMatrix gamma,
//                            Rcpp::List covariates):
//     dispersion_(std::forward<D>(dispersion)),
//     medium_range_dispersion_(std::forward<E>(medium_range_dispersion)),
//     lambda_(lambda),
//     alpha_(alpha),
//     beta_(beta),
//     gamma_(gamma),
//     covariates_(covariates) {}
//
//   template<typename Point, typename Configuration>
//   double compute_papangelou(const Point& point, const Configuration& configuration) const {
//     const auto dalpha(compute_dispersion(dispersion_, point, alpha_.nrow(), configuration));
//     double alpha_dispersion(detail::compute_alpha_dot_dispersion(point, alpha_, dalpha));
//     const auto dgamma(compute_dispersion(medium_range_dispersion_, point, gamma_.nrow(), configuration));
//     double gamma_dispersion(detail::compute_alpha_dot_dispersion(point, gamma_, dgamma));
//     double beta_covariates(detail::compute_beta_dot_covariates(point, beta_, covariates_));
//     return lambda_[get_type(point)] * std::exp(alpha_dispersion + beta_covariates + gamma_dispersion);
//   }
//
// protected:
//   Dispersion dispersion_;
//   MediumRangeDispersion medium_range_dispersion_;
//   Lambda lambda_;
//   Rcpp::NumericMatrix alpha_;
//   Rcpp::NumericMatrix beta_;
//   Rcpp::NumericMatrix gamma_;
//   Im_list_wrapper covariates_;
// };
//
// template<typename Window, typename Dispersion, typename MediumRangeDispersion, typename Lambda>
// class Exponential_family_model_over_window: public Exponential_family_model<Dispersion, MediumRangeDispersion, Lambda> {
// private:
//   using Model = Exponential_family_model<Dispersion, MediumRangeDispersion, Lambda>;
// public:
//   template<typename D, typename E>
//   Exponential_family_model_over_window(const Window& window,
//                                        const Lambda& lambda,
//                                        D&& dispersion,
//                                        E&& medium_range_dispersion,
//                                        Rcpp::NumericMatrix alpha,
//                                        Rcpp::NumericMatrix beta,
//                                        Rcpp::NumericMatrix gamma,
//                                        Rcpp::List covariates):
//     Model(lambda, std::forward<D>(dispersion), std::forward<E>(medium_range_dispersion), alpha, beta, gamma, covariates),
//     window_(window),
//     beta_dot_covariates_maximum_(detail::compute_beta_dot_covariates_maximum(beta, Model::covariates_)),
//     dot_dispersion_maximum_(detail::compute_alpha_dot_dispersion_maximum(alpha, dispersion.get_maximum())) {
//     const auto gamma_dot_dispersion_maximum(detail::compute_alpha_dot_dispersion_maximum(gamma, medium_range_dispersion.get_maximum()));
//     std::transform(dot_dispersion_maximum_.begin(), dot_dispersion_maximum_.end(), gamma_dot_dispersion_maximum.begin(),
//                    dot_dispersion_maximum_.begin(), std::plus<double>());
//   }
//
//   template<typename Point>
//   auto get_log_normalized_bounding_intensity(const Point& point) const {
//     double beta_covariates_maximum(beta_dot_covariates_maximum_[get_type(point)]);
//     double beta_covariates(detail::compute_beta_dot_covariates(point, Model::beta_, Model::covariates_));
//     return beta_covariates - beta_covariates_maximum;
//   }
//
//   // TODO: This should somehow be restricted to window.
//   auto get_integral() const {
//     double integral(0);
//     const auto number_types(Model::beta_.nrow());
//     for(R_xlen_t i(0); i < number_types; ++i) {
//       integral += Model::lambda_[i] * Model::covariates_.get_integral_of_dot(window_, [i, this](double x) {
//         return std::exp(x + dot_dispersion_maximum_[i]);
//       }, Model::beta_(i, Rcpp::_));
//     }
//     return integral;
//   }
//
//   auto sample_point_from_bounding_intensity() const {
//     // Sample type proportionally to the \lambda_i.
//     const auto random_type(Rcpp::sample(Model::lambda_.size(), 1, false, Rcpp::sugar::probs_t(Model::lambda_), false)[0]);
//     while(true) {
//       const auto sample(window_.sample(random_type));
//       if(exp_rand() + get_log_normalized_bounding_intensity(sample) >= 0) {
//         return sample;
//       }
//     }
//   }
//
//   auto get_upper_bound() const {
//     const auto number_types(Model::lambda_.size());
//     std::vector<double> upper_bound(number_types);
//     using size_t = decltype(Model::lambda_.size());
//     for(size_t i(0); i < number_types; ++i) {
//       const auto value(Model::lambda_[i] * std::exp(dot_dispersion_maximum_[i] + beta_dot_covariates_maximum_[i]));
//       if(std::isinf(value)) {
//         Rcpp::stop("Infinite value obtained as the bound to the Papangelou intensity.");
//       }
//       upper_bound[i] = value;
//     }
//     return upper_bound;
//   }
//
//   // TODO: Factorise and make this function more readable.
//   template<typename Point, typename Configuration>
//   void add_to_L_or_U(double exp_mark, const Point& point, Configuration& l, Configuration& l_complement) const {
//     double dot_dispersion_maximum(dot_dispersion_maximum_[get_type(point)]);
//     if(detail::is_alpha_non_negative(point, Model::alpha_) && detail::is_alpha_non_negative(point, Model::gamma_)) {
//       const auto alpha_dispersion_u(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l, l_complement));
//       const auto gamma_dispersion_u(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l, l_complement));
//       const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_u));
//       const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_u));
//       if(alpha_dispersion + gamma_dispersion - dot_dispersion_maximum + exp_mark > 0) {
//         const auto alpha_dispersion_l(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l));
//         const auto gamma_dispersion_l(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l));
//         const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_l));
//         const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_l));
//         add_point(alpha_dispersion + gamma_dispersion - dot_dispersion_maximum + exp_mark > 0 ? l : l_complement, point);
//       }
//     } else if(detail::is_alpha_non_positive(point, Model::alpha_) && detail::is_alpha_non_positive(point, Model::gamma_)) {
//       const auto alpha_dispersion_l(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l));
//       const auto gamma_dispersion_l(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l));
//       const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_l));
//       const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_l));
//       if(alpha_dispersion + gamma_dispersion - dot_dispersion_maximum + exp_mark > 0) {
//         const auto alpha_dispersion_u(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l, l_complement));
//         const auto gamma_dispersion_u(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l, l_complement));
//         const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_u));
//         const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_u));
//         add_point(alpha_dispersion + gamma_dispersion - dot_dispersion_maximum + exp_mark > 0 ? l : l_complement, point);
//       }
//     } else {
//       const auto alpha_dispersion_l(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l));
//       const auto gamma_dispersion_l(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l));
//       const auto alpha_dispersion_u(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l, l_complement));
//       const auto gamma_dispersion_u(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l, l_complement));
//       const auto positive_alpha(detail::compute_positive_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_u));
//       const auto negative_alpha(detail::compute_negative_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_l));
//       const auto positive_gamma(detail::compute_positive_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_u));
//       const auto negative_gamma(detail::compute_negative_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_l));
//       if(positive_alpha + negative_alpha + positive_gamma + negative_gamma - dot_dispersion_maximum + exp_mark > 0) {
//         const auto positive_alpha(detail::compute_positive_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_l));
//         const auto negative_alpha(detail::compute_negative_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_u));
//         const auto positive_gamma(detail::compute_positive_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_l));
//         const auto negative_gamma(detail::compute_negative_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_u));
//         add_point(positive_alpha + negative_alpha + positive_gamma + negative_gamma - dot_dispersion_maximum + exp_mark > 0 ? l : l_complement, point);
//       }
//     }
//   }
//
//   double compute_log_alpha_min_lower_bound(R_xlen_t type) const {
//     double sum_alpha(0);
//     for(R_xlen_t j(0); j < Model::alpha_.ncol(); ++j) {
//       const auto a(Model::alpha_(type, j));
//       if(a > 0) {
//         sum_alpha -= a;
//       } else {
//         sum_alpha += a;
//       }
//     }
//     double sum_gamma(0);
//     for(R_xlen_t j(0); j < Model::gamma_.ncol(); ++j) {
//       const auto g(Model::gamma_(type, j));
//       if(g > 0) {
//         sum_gamma -= g;
//       } else {
//         sum_gamma += g;
//       }
//     }
//
//     return sum_alpha * Model::dispersion_.get_maximum()
//       + sum_gamma * Model::medium_range_dispersion_.get_maximum();
//   }
//
//   const auto& get_window() const {
//     return window_;
//   }
//
// private:
//   Window window_;
//   std::vector<double> beta_dot_covariates_maximum_;
//   std::vector<double> dot_dispersion_maximum_;
// };

template<typename Lambda>
class Truncated_exponential_family_model {
public:
  Truncated_exponential_family_model(const Lambda& lambda,
                                     Rcpp::CharacterVector model,
                                     Rcpp::CharacterVector medium_range_model,
                                     Rcpp::NumericMatrix alpha,
                                     Rcpp::NumericMatrix beta,
                                     Rcpp::NumericMatrix gamma,
                                     Rcpp::List covariates,
                                     unsigned long long int max_points,
                                     Rcpp::NumericMatrix short_range,
                                     Rcpp::NumericMatrix medium_range,
                                     Rcpp::NumericMatrix long_range,
                                     unsigned long long int saturation):
    dispersion_(Saturated_model(model, short_range, saturation)),
    medium_range_dispersion_(Saturated_model(medium_range_model, medium_range, long_range, saturation)),
    lambda_(lambda),
    alpha_(alpha),
    beta_(beta),
    gamma_(gamma),
    covariates_(covariates),
    max_points_(max_points) {}

  template<typename Point, typename Configuration>
  double compute_papangelou(const Point& point, const Configuration& configuration) const {
    if(static_cast<unsigned long long int>(size(configuration)) + 1u <= max_points_) {
      const auto dalpha(compute_dispersion(dispersion_, point, alpha_.nrow(), configuration));
      double alpha_dispersion(detail::compute_alpha_dot_dispersion(point, alpha_, dalpha));
      const auto dgamma(compute_dispersion(medium_range_dispersion_, point, gamma_.nrow(), configuration));
      double gamma_dispersion(detail::compute_alpha_dot_dispersion(point, gamma_, dgamma));
      double beta_covariates(detail::compute_beta_dot_covariates(point, beta_, covariates_));
      return lambda_[get_type(point)] * std::exp(beta_covariates + alpha_dispersion + gamma_dispersion);
    } else {
      return 0.;
    }
  }

protected:
  Saturated_model dispersion_;
  Saturated_model medium_range_dispersion_;
  Lambda lambda_;
  Rcpp::NumericMatrix alpha_;
  Rcpp::NumericMatrix beta_;
  Rcpp::NumericMatrix gamma_;
  Im_list_wrapper covariates_;
  unsigned long long int max_points_;
};

template<typename Lambda>
class Truncated_exponential_family_model_over_window: public Truncated_exponential_family_model<Lambda> {
private:
  using Model = Truncated_exponential_family_model<Lambda>;

  template<typename D>
  auto get_dispersion_maximum(const D& dispersion) const {
    return dispersion.get_maximum() + static_cast<double>(Model::max_points_);
  }
public:
  Truncated_exponential_family_model_over_window(const Window& window,
                                                 const Lambda& lambda,
                                                 Rcpp::CharacterVector model,
                                                 Rcpp::CharacterVector medium_range_model,
                                                 Rcpp::NumericMatrix alpha,
                                                 Rcpp::NumericMatrix beta,
                                                 Rcpp::NumericMatrix gamma,
                                                 Rcpp::List covariates,
                                                 unsigned long long int max_points,
                                                 Rcpp::NumericMatrix short_range,
                                                 Rcpp::NumericMatrix medium_range,
                                                 Rcpp::NumericMatrix long_range,
                                                 unsigned long long int saturation):
    Model(lambda, model, medium_range_model, alpha, beta, gamma, covariates, max_points, short_range, medium_range, long_range, saturation),
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

  // TODO: This should somehow be restricted to window.
  auto get_integral() const {
    double integral(0);
    const auto number_types(Model::beta_.nrow());
    for(R_xlen_t i(0); i < number_types; ++i) {
      integral += Model::lambda_[i] * Model::covariates_.get_integral_of_dot(window_, [i, this](double x) {
        return std::exp(x + dot_dispersion_maximum_[i]);
      }, Model::beta_(i, Rcpp::_));
    }
    return integral;
  }

  auto sample_point_from_bounding_intensity() const {
    // Sample type proportionally to the \  lambda_i.
    const auto random_type(Rcpp::sample(Model::lambda_.size(), 1, false, Rcpp::sugar::probs_t(Model::lambda_), false)[0]);
    while(true) {
      const auto sample(window_.sample(random_type));
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
      const auto value(Model::lambda_[i] * std::exp(dot_dispersion_maximum_[i] + beta_dot_covariates_maximum_[i]));
      if(std::isinf(value)) {
        Rcpp::stop("Infinite value obtained as the bound to the Papangelou intensity.");
      }
      upper_bound[i] = value;
    }
    return upper_bound;
  }

  // TODO: Factorise and make this function more readable.
  // TODO: Doesn't work for truncated version--fix.
  // TODO: The IPPP in the CFTP algorithm should be conditioned to have <= saturation points.
  template<typename Point, typename Configuration>
  void add_to_L_or_U(double exp_mark, const Point& point, Configuration& l, Configuration& l_complement) const {
    // TODO: Continue thinking about whether or not this check is correct. Write properly what this function is doing.
    if(size(l) + size(l_complement) >= Model::max_points_) {
      return;
    }
    double dot_dispersion_maximum(dot_dispersion_maximum_[get_type(point)]);
    if(detail::is_alpha_non_negative(point, Model::alpha_) && detail::is_alpha_non_negative(point, Model::gamma_)) {
      const auto alpha_dispersion_u(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l, l_complement));
      const auto gamma_dispersion_u(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l, l_complement));
      const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_u));
      const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_u));
      if(alpha_dispersion + gamma_dispersion - dot_dispersion_maximum + exp_mark > 0) {
        const auto alpha_dispersion_l(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l));
        const auto gamma_dispersion_l(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l));
        const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_l));
        const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_l));
          add_point(alpha_dispersion + gamma_dispersion - dot_dispersion_maximum + exp_mark > 0 ? l : l_complement, point);
      }
    } else if(detail::is_alpha_non_positive(point, Model::alpha_) && detail::is_alpha_non_positive(point, Model::gamma_)) {
      const auto alpha_dispersion_l(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l));
      const auto gamma_dispersion_l(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l));
      const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_l));
      const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_l));
      if(alpha_dispersion + gamma_dispersion - dot_dispersion_maximum + exp_mark > 0) {
        const auto alpha_dispersion_u(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l, l_complement));
        const auto gamma_dispersion_u(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l, l_complement));
        const auto alpha_dispersion(detail::compute_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_u));
        const auto gamma_dispersion(detail::compute_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_u));
        add_point(alpha_dispersion + gamma_dispersion - dot_dispersion_maximum + exp_mark > 0 ? l : l_complement, point);
      }
    } else {
      const auto alpha_dispersion_l(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l));
      const auto gamma_dispersion_l(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l));
      const auto alpha_dispersion_u(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l, l_complement));
      const auto gamma_dispersion_u(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l, l_complement));
      const auto positive_alpha(detail::compute_positive_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_u));
      const auto negative_alpha(detail::compute_negative_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_l));
      const auto positive_gamma(detail::compute_positive_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_u));
      const auto negative_gamma(detail::compute_negative_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_l));
      if(positive_alpha + negative_alpha + positive_gamma + negative_gamma - dot_dispersion_maximum + exp_mark > 0) {
        const auto positive_alpha(detail::compute_positive_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_l));
        const auto negative_alpha(detail::compute_negative_alpha_dot_dispersion(point, Model::alpha_, alpha_dispersion_u));
        const auto positive_gamma(detail::compute_positive_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_l));
        const auto negative_gamma(detail::compute_negative_alpha_dot_dispersion(point, Model::gamma_, gamma_dispersion_u));
        add_point(positive_alpha + negative_alpha + positive_gamma + negative_gamma - dot_dispersion_maximum + exp_mark > 0 ? l : l_complement, point);
      }
    }
  }

  // TODO: Obviously don't want to use this anymore since fixed version is not lower-bounded.
  double compute_log_alpha_min_lower_bound(R_xlen_t) const {
    return -std::numeric_limits<double>::infinity();
    // double sum_alpha(0);
    // for(R_xlen_t j(0); j < Model::alpha_.ncol(); ++j) {
    //   const auto a(Model::alpha_(type, j));
    //   if(a > 0) {
    //     sum_alpha -= a;
    //   } else {
    //     sum_alpha += a;
    //   }
    // }
    // double sum_gamma(0);
    // for(R_xlen_t j(0); j < Model::gamma_.ncol(); ++j) {
    //   const auto g(Model::gamma_(type, j));
    //   if(g > 0) {
    //     sum_gamma -= g;
    //   } else {
    //     sum_gamma += g;
    //   }
    // }
    //
    // return sum_alpha * Model::dispersion_.get_maximum()
    //   + sum_gamma * Model::medium_range_dispersion_.get_maximum();
  }

  const auto& get_window() const {
    return window_;
  }

private:
  Window window_;
  std::vector<double> beta_dot_covariates_maximum_;
  std::vector<double> dot_dispersion_maximum_;
};

template<typename F, typename Lambda, typename... Args>
inline auto call_on_model(Rcpp::CharacterVector model,
                          Rcpp::CharacterVector medium_range_model,
                          const Lambda& lambda,
                          Rcpp::NumericMatrix short_range,
                          Rcpp::NumericMatrix medium_range,
                          Rcpp::NumericMatrix long_range,
                          unsigned long long int saturation,
                          unsigned long long int max_points,
                          const F& f,
                          Args... args) {
  const auto dispersion(Saturated_model(model, short_range, saturation));
  const auto medium_range_dispersion(Saturated_model(medium_range_model, medium_range, long_range, saturation));
  using Model_type = Truncated_exponential_family_model<Lambda>;
  return f(Model_type(lambda, std::move(dispersion), std::move(medium_range_dispersion), args..., max_points));
}

template<typename F, typename Lambda, typename... Args>
inline auto call_on_model(const Window& window,
                          Rcpp::CharacterVector model,
                          Rcpp::CharacterVector medium_range_model,
                          const Lambda& lambda,
                          Rcpp::NumericMatrix short_range,
                          Rcpp::NumericMatrix medium_range,
                          Rcpp::NumericMatrix long_range,
                          unsigned long long int saturation,
                          unsigned long long int max_points,
                          const F& f,
                          Args... args) {
  const auto dispersion(Saturated_model(model, short_range, saturation));
  const auto medium_range_dispersion(Saturated_model(medium_range_model, medium_range, long_range, saturation));
  using Model_type = Truncated_exponential_family_model_over_window<Lambda>;
  // TODO: Fix the endpoints
  return f(Model_type(window, lambda, std::move(dispersion), std::move(medium_range_dispersion), args..., max_points));
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_EXPONENTIAL_FAMILY
