#ifndef INCLUDE_PPJSDM_EXPONENTIAL_FAMILY
#define INCLUDE_PPJSDM_EXPONENTIAL_FAMILY

#include <Rcpp.h>

#include "saturated_model.hpp"
#include "../configuration/configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"
#include "../simulation/inhomogeneous_ppp.hpp"
#include "../utility/algebra.hpp"
#include "../utility/approximate_draw.hpp"
#include "../utility/im_wrapper.hpp"

#include <algorithm> // std::transform
#include <cmath> // std::exp, std::isinf
#include <functional> // std::plus
#include <vector> // std::vector

namespace ppjsdm {
namespace detail {

template<typename Point, typename Beta, typename Covariates>
inline auto compute_beta_dot_covariates(const Point& point, const Beta& beta, const Covariates& covariates) {
  double sum(0);
  const auto number_covariates(covariates.size());
  using size_t = std::remove_cv_t<decltype(covariates.size())>;
  for(size_t i(0); i < number_covariates; ++i) {
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
  const auto number_covariates(beta.nrow());
  std::vector<double> result(number_covariates);
  if(covariates.size() > 0) {
    using size_t = std::remove_cv_t<decltype(number_covariates)>;
    for(size_t i(0); i < number_covariates; ++i) {
      result[i] = compute_beta_dot_covariates_maximum(i, beta, covariates);
    }
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
  auto compute_log_papangelou(const Point& point, const Configuration& configuration) const {
    auto inner_product(beta0_[get_type(point)]);
    if(!is_column_zero(alpha_, get_type(point))) {
      const auto dispersion(compute_dispersion(dispersion_, point, alpha_.nrow(), configuration));
      inner_product += matrix_times_vector_at_index(alpha_, dispersion, get_type(point));
    }
    if(!is_column_zero(gamma_, get_type(point))) {
      const auto dispersion(compute_dispersion(medium_range_dispersion_, point, gamma_.nrow(), configuration));
      inner_product += matrix_times_vector_at_index(gamma_, dispersion, get_type(point));
    }
    inner_product += detail::compute_beta_dot_covariates(point, beta_, covariates_);
    return inner_product;
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
    dot_dispersion_maximum_(positive_matrix_times_vector(alpha, get_dispersion_maximum(Model::dispersion_))) {
    const auto gamma_dot_dispersion_maximum(positive_matrix_times_vector(gamma, get_dispersion_maximum(Model::medium_range_dispersion_)));
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
    const auto number_types(Model::get_number_types());
    using size_t = std::remove_cv_t<decltype(number_types)>;
    for(size_t type(0); type < number_types; ++type) {
      integral += Model::covariates_.get_integral_of_dot(window_, [type, this](double x) {
        return std::exp(x + Model::beta0_[type] + dot_dispersion_maximum_[type]);
      }, Model::beta_(type, Rcpp::_));
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
    const auto number_types(Model::get_number_types());
    std::vector<double> upper_bound(number_types);
    using size_t = std::remove_cv_t<decltype(number_types)>;
    for(size_t type(0); type < number_types; ++type) {
      const auto value(Model::beta0_[type] + beta_dot_covariates_maximum_[type]);
      if(std::isinf(value)) {
        Rcpp::Rcout << Model::beta0_[type] << '\n';
        Rcpp::Rcout << beta_dot_covariates_maximum_[type] << '\n';
        Rcpp::stop("Infinite value obtained as the bound to the approximate PPP intensity.");
      }
      upper_bound[type] = std::exp(value);
    }
    return upper_bound;
  }

  auto get_upper_bound() const {
    const auto number_types(Model::get_number_types());
    std::vector<double> upper_bound(number_types);
    using size_t = std::remove_cv_t<decltype(number_types)>;
    for(size_t type(0); type < number_types; ++type) {
      const auto value(Model::beta0_[type] + dot_dispersion_maximum_[type] + beta_dot_covariates_maximum_[type]);
      if(std::isinf(value)) {
        Rcpp::Rcout << Model::beta0_[type] << '\n';
        Rcpp::Rcout << dot_dispersion_maximum_[type] << '\n';
        Rcpp::Rcout << beta_dot_covariates_maximum_[type] << '\n';
        Rcpp::stop("Infinite value obtained as the bound to the Papangelou intensity.");
      }
      upper_bound[type] = std::exp(value);
    }
    return upper_bound;
  }

  // TODO: Factorise and make this function more readable.
  template<typename Point, typename Configuration>
  void add_to_L_or_U(double exp_mark, const Point& point, Configuration& l, Configuration& l_complement) const {
    // TODO: Write properly what this function is doing.
    double dot_dispersion_maximum(dot_dispersion_maximum_[get_type(point)]);
    if(is_column_nonnegative(Model::alpha_, get_type(point)) && is_column_nonnegative(Model::gamma_, get_type(point))) {
      auto log_alpha(-dot_dispersion_maximum);
      if(!is_column_zero(Model::alpha_, get_type(point))) {
        const auto short_range_dispersion(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l, l_complement));
        log_alpha += matrix_times_vector_at_index(Model::alpha_, short_range_dispersion, get_type(point));
      }
      if(!is_column_zero(Model::gamma_, get_type(point))) {
        const auto medium_range_dispersion(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l, l_complement));
        log_alpha += matrix_times_vector_at_index(Model::gamma_, medium_range_dispersion, get_type(point));
      }
      if(log_alpha > 0.) {
        Rcpp::stop("Bad upper-bound in the computation of alpha in the CFTP algorithm.");
      }
      if(log_alpha + exp_mark > 0.) {
        auto log_alpha(-dot_dispersion_maximum);
        if(!is_column_zero(Model::alpha_, get_type(point))) {
          const auto short_range_dispersion(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l));
          log_alpha += matrix_times_vector_at_index(Model::alpha_, short_range_dispersion, get_type(point));
        }
        if(!is_column_zero(Model::gamma_, get_type(point))) {
          const auto medium_range_dispersion(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l));
          log_alpha += matrix_times_vector_at_index(Model::gamma_, medium_range_dispersion, get_type(point));
        }
        if(log_alpha > 0.) {
          Rcpp::stop("Bad upper-bound in the computation of alpha in the CFTP algorithm.");
        }
        add_point(log_alpha + exp_mark > 0. ? l : l_complement, point);
      }
    } else if(is_column_nonpositive(Model::alpha_, get_type(point)) && is_column_nonpositive(Model::gamma_, get_type(point))) {
      auto log_alpha(-dot_dispersion_maximum);
      if(!is_column_zero(Model::alpha_, get_type(point))) {
        const auto short_range_dispersion(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l));
        log_alpha += matrix_times_vector_at_index(Model::alpha_, short_range_dispersion, get_type(point));
      }
      if(!is_column_zero(Model::gamma_, get_type(point))) {
        const auto medium_range_dispersion(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l));
        log_alpha += matrix_times_vector_at_index(Model::gamma_, medium_range_dispersion, get_type(point));
      }
      if(log_alpha > 0.) {
        Rcpp::stop("Bad upper-bound in the computation of alpha in the CFTP algorithm.");
      }
      if(log_alpha + exp_mark > 0.) {
        auto log_alpha(-dot_dispersion_maximum);
        if(!is_column_zero(Model::alpha_, get_type(point))) {
          const auto short_range_dispersion(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l, l_complement));
          log_alpha += matrix_times_vector_at_index(Model::alpha_, short_range_dispersion, get_type(point));
        }
        if(!is_column_zero(Model::gamma_, get_type(point))) {
          const auto medium_range_dispersion(compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l, l_complement));
          log_alpha += matrix_times_vector_at_index(Model::gamma_, medium_range_dispersion, get_type(point));
        }
        if(log_alpha > 0.) {
          Rcpp::stop("Bad upper-bound in the computation of alpha in the CFTP algorithm.");
        }
        add_point(log_alpha + exp_mark > 0. ? l : l_complement, point);
      }
    } else {
      using dispersion_t = std::remove_cv_t<std::remove_reference_t<decltype(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l))>>;
      dispersion_t short_range_dispersion_l, short_range_dispersion_u, medium_range_dispersion_l, medium_range_dispersion_u;
      auto log_alpha(-dot_dispersion_maximum);
      if(!is_column_zero(Model::alpha_, get_type(point))) {
        short_range_dispersion_u = compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l, l_complement);
        short_range_dispersion_l = compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l);
        log_alpha += positive_matrix_times_vector_at_index(Model::alpha_, short_range_dispersion_u, get_type(point));
        log_alpha += negative_matrix_times_vector_at_index(Model::alpha_, short_range_dispersion_l, get_type(point));
      }
      if(!is_column_zero(Model::gamma_, get_type(point))) {
        medium_range_dispersion_u = compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l, l_complement);
        medium_range_dispersion_l = compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l);
        log_alpha += positive_matrix_times_vector_at_index(Model::gamma_, medium_range_dispersion_u, get_type(point));
        log_alpha += negative_matrix_times_vector_at_index(Model::gamma_, medium_range_dispersion_l, get_type(point));
      }
      if(log_alpha > 0.) {
        Rcpp::stop("Bad upper-bound in the computation of alpha in the CFTP algorithm.");
      }
      if(log_alpha + exp_mark > 0.) {
        auto log_alpha(-dot_dispersion_maximum);
        if(!is_column_zero(Model::alpha_, get_type(point))) {
          log_alpha += positive_matrix_times_vector_at_index(Model::alpha_, short_range_dispersion_l, get_type(point));
          log_alpha += negative_matrix_times_vector_at_index(Model::alpha_, short_range_dispersion_u, get_type(point));
        }
        if(!is_column_zero(Model::gamma_, get_type(point))) {
          log_alpha += positive_matrix_times_vector_at_index(Model::gamma_, medium_range_dispersion_l, get_type(point));
          log_alpha += negative_matrix_times_vector_at_index(Model::gamma_, medium_range_dispersion_u, get_type(point));
        }
        if(log_alpha > 0.) {
          Rcpp::stop("Bad upper-bound in the computation of alpha in the CFTP algorithm.");
        }
        add_point(log_alpha + exp_mark > 0. ? l : l_complement, point);
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
                                                                         return model.get_log_normalized_bounding_intensity(point);
                                                                       },
                                                                       model.get_upper_bound_approximate_ppp_intensity(),
                                                                       model.get_number_types()));
    return configuration;
  }
};

} // namespace detail

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_EXPONENTIAL_FAMILY
