#ifndef INCLUDE_PPJSDM_EXPONENTIAL_FAMILY
#define INCLUDE_PPJSDM_EXPONENTIAL_FAMILY

#include <Rcpp.h>

#include "compute_dispersion.hpp"
#include "compute_dispersion_fitting.hpp"
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
    return beta0_[get_type(point)] +
      compute_total_dispersion(point, configuration) +
      detail::compute_beta_dot_covariates(point, beta_, covariates_);
  }

  template<typename Point, typename Configuration>
  double compute_papangelou(const Point& point, const Configuration& configuration) const {
    return std::exp(compute_log_papangelou(point, configuration));
  }

  // TODO: Code quality is awful for this, make a unique compute_papangelou function, same for compute_total_dispersion
  template<typename Points, typename Configuration>
  auto compute_papangelou_vectorized(const Points& points, const Configuration& configuration) const {
    auto dispersion(compute_total_dispersion_vectorized(points, configuration));
    for(typename decltype(dispersion)::size_type i(0); i < dispersion.size(); ++i) {
      dispersion[i] += beta0_[get_type(points[i])];
      dispersion[i] += detail::compute_beta_dot_covariates(points[i], beta_, covariates_);
      dispersion[i] = std::exp(dispersion[i]);
    }
    return dispersion;
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

  template<typename Point, typename... Configurations>
  auto compute_total_dispersion(const Point& point, Configurations&... configurations) const {
    using return_t = decltype(matrix_times_vector_at_index(alpha_, compute_dispersion(dispersion_, point, alpha_.nrow(), configurations...), 0));
    const auto number_types(alpha_.nrow());
    return_t return_value(0.);
    if(!is_column_zero(alpha_, get_type(point))) {
      const auto short_range_dispersion(compute_dispersion(dispersion_, point, number_types, configurations...));
      return_value += matrix_times_vector_at_index(alpha_, short_range_dispersion, get_type(point));
    }
    if(!is_column_zero(gamma_, get_type(point))) {
      const auto medium_range_dispersion(compute_dispersion(medium_range_dispersion_, point, number_types, configurations...));
      return_value += matrix_times_vector_at_index(gamma_, medium_range_dispersion, get_type(point));
    }
    return return_value;
  }

  template<typename Points, typename... Configurations>
  auto compute_total_dispersion_vectorized(const Points& points, Configurations&... configurations) const {
    const auto number_types(alpha_.nrow());
    std::vector<decltype(matrix_times_vector_at_index(alpha_, compute_dispersion_for_fitting<false>(dispersion_, number_types, configurations..., points)[0], 0))> return_value(points.size());
    if(!is_column_zero(alpha_, get_type(points[0]))) {
      const auto short_range_dispersion(compute_dispersion_for_fitting<false>(dispersion_, number_types, configurations..., points));
      for(typename Points::size_type i(0); i < return_value.size(); ++i) {
        return_value[i] += matrix_times_vector_at_index(alpha_, short_range_dispersion[i], get_type(points[i]));
      }
    }
    if(!is_column_zero(gamma_, get_type(points[0]))) {
      const auto medium_range_dispersion(compute_dispersion_for_fitting<false>(medium_range_dispersion_, number_types, configurations..., points));
      for(typename Points::size_type i(0); i < return_value.size(); ++i) {
        return_value[i] += matrix_times_vector_at_index(gamma_, medium_range_dispersion[i], get_type(points[i]));
      }
    }
    return return_value;
  }
};

template<typename Lambda>
class Truncated_exponential_family_model_over_window: public Truncated_exponential_family_model<Lambda> {
private:
  using Model = Truncated_exponential_family_model<Lambda>;

  template<typename D>
  auto get_dispersion_maximum(const D& dispersion) const {
    if(dispersion.is_nonincreasing_after_lower_endpoint()) {
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
      const auto value(std::exp(Model::beta0_[type] + beta_dot_covariates_maximum_[type]));
      if(std::isinf(value)) {
        Rcpp::Rcout << Model::beta0_[type] << '\n';
        Rcpp::Rcout << beta_dot_covariates_maximum_[type] << '\n';
        Rcpp::stop("Infinite value obtained as the bound to the approximate PPP intensity.");
      }
      upper_bound[type] = value;
    }
    return upper_bound;
  }

  auto get_upper_bound() const {
    const auto number_types(Model::get_number_types());
    std::vector<double> upper_bound(number_types);
    using size_t = std::remove_cv_t<decltype(number_types)>;
    for(size_t type(0); type < number_types; ++type) {
      const auto value(std::exp(Model::beta0_[type] + dot_dispersion_maximum_[type] + beta_dot_covariates_maximum_[type]));
      if(std::isinf(value)) {
        Rcpp::Rcout << Model::beta0_[type] << '\n';
        Rcpp::Rcout << dot_dispersion_maximum_[type] << '\n';
        Rcpp::Rcout << beta_dot_covariates_maximum_[type] << '\n';
        Rcpp::stop("Infinite value obtained as the bound to the Papangelou intensity.");
      }
      upper_bound[type] = value;
    }
    return upper_bound;
  }

  // TODO: Factorise and make this function more readable.
  // TODO: Write properly what this function is doing.
  template<typename Point, typename Configuration>
  void add_to_L_or_U(double exp_mark, const Point& point, Configuration& l, Configuration& l_complement) const {
    const auto dot_dispersion_maximum(dot_dispersion_maximum_[get_type(point)]);
    using FloatType = std::remove_cv_t<decltype(dot_dispersion_maximum)>;
    FloatType alpha_need_to_add(0.), alpha_which_one_to_add_to(0.);
    if(is_column_nonnegative(Model::alpha_, get_type(point)) && is_column_nonnegative(Model::gamma_, get_type(point))) {
      alpha_need_to_add = Model::compute_total_dispersion(point, l, l_complement) - dot_dispersion_maximum;
      if(alpha_need_to_add + exp_mark > 0.) {
        alpha_which_one_to_add_to = Model::compute_total_dispersion(point, l) - dot_dispersion_maximum;
        add_point(alpha_which_one_to_add_to + exp_mark > 0. ? l : l_complement, point);
      }
    } else if(is_column_nonpositive(Model::alpha_, get_type(point)) && is_column_nonpositive(Model::gamma_, get_type(point))) {
      alpha_need_to_add = Model::compute_total_dispersion(point, l) - dot_dispersion_maximum;
      if(alpha_need_to_add + exp_mark > 0.) {
        alpha_which_one_to_add_to = Model::compute_total_dispersion(point, l, l_complement) - dot_dispersion_maximum;
        add_point(alpha_which_one_to_add_to + exp_mark > 0. ? l : l_complement, point);
      }
    } else {
      using dispersion_t = std::remove_cv_t<std::remove_reference_t<decltype(compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l))>>;
      dispersion_t short_range_dispersion_l, short_range_dispersion_u, medium_range_dispersion_l, medium_range_dispersion_u;
      alpha_need_to_add = -dot_dispersion_maximum;
      if(!is_column_zero(Model::alpha_, get_type(point))) {
        // TODO: These two can definitely be computed faster by doing both at once.
        short_range_dispersion_u = compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l, l_complement);
        short_range_dispersion_l = compute_dispersion(Model::dispersion_, point, Model::alpha_.nrow(), l);
        alpha_need_to_add += positive_matrix_times_vector_at_index(Model::alpha_, short_range_dispersion_u, get_type(point));
        alpha_need_to_add += negative_matrix_times_vector_at_index(Model::alpha_, short_range_dispersion_l, get_type(point));
      }
      if(!is_column_zero(Model::gamma_, get_type(point))) {
        // TODO: These two can definitely be computed faster by doing both at once.
        medium_range_dispersion_u = compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l, l_complement);
        medium_range_dispersion_l = compute_dispersion(Model::medium_range_dispersion_, point, Model::gamma_.nrow(), l);
        alpha_need_to_add += positive_matrix_times_vector_at_index(Model::gamma_, medium_range_dispersion_u, get_type(point));
        alpha_need_to_add += negative_matrix_times_vector_at_index(Model::gamma_, medium_range_dispersion_l, get_type(point));
      }
      if(alpha_need_to_add + exp_mark > 0.) {
        alpha_which_one_to_add_to = -dot_dispersion_maximum;
        if(!is_column_zero(Model::alpha_, get_type(point))) {
          alpha_which_one_to_add_to += positive_matrix_times_vector_at_index(Model::alpha_, short_range_dispersion_l, get_type(point));
          alpha_which_one_to_add_to += negative_matrix_times_vector_at_index(Model::alpha_, short_range_dispersion_u, get_type(point));
        }
        if(!is_column_zero(Model::gamma_, get_type(point))) {
          alpha_which_one_to_add_to += positive_matrix_times_vector_at_index(Model::gamma_, medium_range_dispersion_l, get_type(point));
          alpha_which_one_to_add_to += negative_matrix_times_vector_at_index(Model::gamma_, medium_range_dispersion_u, get_type(point));
        }
        add_point(alpha_which_one_to_add_to + exp_mark > 0. ? l : l_complement, point);
      }
    }
    if(alpha_need_to_add > 0. || alpha_which_one_to_add_to > 0.) {
      Rcpp::stop("Internal error: Incorrect alpha value, probably due to bad upper-bound to dispersion in the CFTP algorithm.");
    }
  }

  double compute_log_alpha_min_lower_bound(R_xlen_t type) const {
    double sum_alpha(0.);
    double sum_gamma(0.);
    for(R_xlen_t other_type(0); other_type < Model::alpha_.ncol(); ++other_type) {
      sum_alpha -= std::abs(Model::alpha_(type, other_type));
      sum_gamma -= std::abs(Model::gamma_(type, other_type));
    }
    return sum_alpha  * get_dispersion_maximum(Model::dispersion_) +
      sum_gamma * get_dispersion_maximum(Model::medium_range_dispersion_);
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
