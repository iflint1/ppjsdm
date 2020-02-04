#ifndef INCLUDE_PPJSDM_EXPONENTIAL_FAMILY
#define INCLUDE_PPJSDM_EXPONENTIAL_FAMILY

#include <Rcpp.h>

#include "saturated_varphi_model.hpp"
#include "../point/point_manipulation.hpp"
#include "../utility/im_wrapper.hpp"

#include <cmath> // std::exp, std::log, std::fabs, std::isinf
#include <type_traits> // std::remove_const, std::remove_reference, std::is_same, std::enable_if
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

template<typename Dispersion, typename Lambda, typename Alpha, typename Beta>
class Exponential_family_model {
public:
  template<typename D, std::enable_if_t<std::is_same<Dispersion, std::remove_reference_t<D>>::value>* = nullptr>
  Exponential_family_model(const Lambda& lambda,
                           const Alpha& alpha,
                           const Beta& beta,
                           Rcpp::List covariates,
                           D&& dispersion):
    dispersion_(std::forward<D>(dispersion)),
    lambda_(lambda),
    alpha_(alpha),
    beta_(beta),
    covariates_(covariates) {}

  class Normalised_dominating_intensity {
  public:
    template<typename V>
    Normalised_dominating_intensity(const Lambda& lambda,
                                    const Beta& beta,
                                    const Im_list_wrapper& covariates,
                                    V&& alpha_dot_dispersion_maximum):
    lambda_(lambda),
    beta_(beta),
    covariates_(covariates),
    beta_dot_covariates_maximum_(detail::compute_beta_dot_covariates_maximum(beta, covariates)),
    alpha_dot_dispersion_maximum_(std::forward<V>(alpha_dot_dispersion_maximum)) {}

    template<typename Point>
    auto operator()(const Point& point) const {
      double beta_covariates_maximum(beta_dot_covariates_maximum_[get_type(point)]);
      double beta_covariates(detail::compute_beta_dot_covariates(point, beta_, covariates_));
      return std::exp(beta_covariates - beta_covariates_maximum);
    }

    auto get_integral() const {
      double integral(0);
      const auto number_types(beta_.nrow());
      for(R_xlen_t i(0); i < number_types; ++i) {
        integral += covariates_.get_integral_of_dot([](double x) { return std::exp(x); }, beta_(i, Rcpp::_));
      }
      return integral / static_cast<double>(number_types);
    }

    auto get_upper_bound() const {
      const auto number_types(lambda_.size());
      std::vector<double> upper_bound(number_types);
      using size_t = decltype(lambda_.size());
      for(size_t i(0); i < number_types; ++i) {
        const auto value(lambda_[i] * std::exp(alpha_dot_dispersion_maximum_[i] + beta_dot_covariates_maximum_[i]));
        if(std::isinf(value)) {
          Rcpp::stop("Infinite value obtained as the bound to the Papangelou intensity.");
        }
        upper_bound[i] = value;
      }
      return upper_bound;
    }
  private:
    Lambda lambda_;
    Beta beta_;
    Im_list_wrapper covariates_;
    std::vector<double> beta_dot_covariates_maximum_;
    std::vector<double> alpha_dot_dispersion_maximum_;
  };

  template<typename Window>
  auto get_normalised_dominating_intensity(const Window& window) const {
    return Normalised_dominating_intensity(lambda_, beta_, covariates_, detail::compute_alpha_dot_dispersion_maximum(alpha_, dispersion_.get_maximum(window)));
  }

  template<typename Point, typename... Configurations>
  double compute_papangelou(const Point& point, Configurations&&... configurations) const {
    double alpha_dispersion(detail::compute_alpha_dot_dispersion(point, alpha_, dispersion_, std::forward<Configurations>(configurations)...));
    double beta_covariates(detail::compute_beta_dot_covariates(point, beta_, covariates_));
    return lambda_[get_type(point)] * std::exp(alpha_dispersion + beta_covariates);
  }

  template<typename Window, typename Point, typename... Configurations>
  double compute_normalised_papangelou(const Window& window, const Point& point, Configurations&&... configurations) const {
    double alpha_dispersion_maximum(detail::compute_alpha_dot_dispersion_maximum(alpha_, dispersion_.get_maximum(window), get_type(point)));
    double alpha_dispersion(detail::compute_alpha_dot_dispersion(point, alpha_, dispersion_, std::forward<Configurations>(configurations)...));
    return std::exp(alpha_dispersion - alpha_dispersion_maximum);
  }

private:
  Dispersion dispersion_;
  Lambda lambda_;
  Alpha alpha_;
  Beta beta_;
  Im_list_wrapper covariates_;
};

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
    return f(Model_type(lambda, alpha, beta, covariates, std::forward<decltype(varphi)>(varphi)));
  });
}


} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_EXPONENTIAL_FAMILY
