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

#include <algorithm> // std::all_of, std::fill, std::transform
#include <cmath> // std::exp, std::isinf
#include <functional> // std::plus
#include <random> // Distributions
#include <string> // std::string
#include <vector> // std::vector

namespace ppjsdm {
namespace detail {

inline auto convert_list_of_numeric_matrices_to_vector(Rcpp::List numeric_matrices) {
  using conversion_t = ppjsdm::Lightweight_matrix<double>;
  std::vector<conversion_t> converted_matrices(numeric_matrices.size());
  for(decltype(converted_matrices.size()) i(0); i < converted_matrices.size(); ++i) {
    converted_matrices[i] = conversion_t(Rcpp::as<Rcpp::NumericMatrix>(numeric_matrices[i]));
  }
  return converted_matrices;
}

template<typename Point, typename Beta, typename Covariates>
inline auto compute_beta_dot_covariates(const Point& point, const Beta& beta, const Covariates& covariates) {
  decltype(beta(0, 0) * covariates[0](point)) sum(0);
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
  std::vector<decltype(compute_beta_dot_covariates_maximum(0, beta, covariates))> result(number_covariates);
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
                                     Rcpp::List model,
                                     Rcpp::CharacterVector medium_range_model,
                                     Rcpp::List alpha,
                                     Rcpp::NumericMatrix beta,
                                     Rcpp::NumericMatrix gamma,
                                     Rcpp::List covariates,
                                     Rcpp::List short_range,
                                     Rcpp::NumericMatrix medium_range,
                                     Rcpp::NumericMatrix long_range,
                                     unsigned long long int saturation):
    medium_range_dispersion_(medium_range_model, medium_range, long_range, saturation),
    beta0_(beta0),
    alpha_(detail::convert_list_of_numeric_matrices_to_vector(alpha)),
    beta_(beta),
    gamma_(gamma),
    covariates_(covariates) {
    if(alpha.size() != short_range.size()) {
      Rcpp::stop("Supplied short_range and alpha lists do not have compatible sizes.");
    } else {
      // TODO: Avoid emplace_back
      dispersion_ = std::vector<Saturated_model<double>>{};
      for(decltype(alpha.size()) i(0); i < alpha.size(); ++i) {
        dispersion_.emplace_back(Saturated_model<double>(model[i], short_range[i], saturation));
      }
    }
  }

  // TODO: Code quality is awful for this, make a unique compute_papangelou function, same for compute_total_dispersion
  template<typename Points, typename Configuration>
  auto compute_papangelou_vectorized(const Points& points, const Configuration& configuration, int nthreads) const {
    auto dispersion(compute_total_dispersion_vectorized(points, configuration, nthreads));
    for(typename decltype(dispersion)::size_type i(0); i < dispersion.size(); ++i) {
      dispersion[i] += beta0_[get_type(points[i])];
      dispersion[i] += detail::compute_beta_dot_covariates(points[i], beta_, covariates_);
      dispersion[i] = std::exp(dispersion[i]);
    }
    return dispersion;
  }

  template<typename Point, typename... Configurations>
  auto compute_papangelou(const Point& point, const Configurations&... configurations) const {
    return std::exp(beta0_[get_type(point)] +
                    compute_total_dispersion(point, configurations...) +
                    detail::compute_beta_dot_covariates(point, beta_, covariates_));
  }

  auto get_number_types() const {
    return beta0_.size();
  }

  auto get_saturation() const {
    return dispersion_[0].get_saturation();
  }

protected:
  std::vector<Saturated_model<double>> dispersion_;
  Saturated_model<double> medium_range_dispersion_;
  Lambda beta0_;
  // TODO: use Lightweight_square_matrix
  std::vector<ppjsdm::Lightweight_matrix<double>> alpha_;
  ppjsdm::Lightweight_matrix<double> beta_;
  ppjsdm::Lightweight_matrix<double> gamma_;
  Im_list_wrapper covariates_;

  template<typename Point, typename... Configurations>
  auto compute_total_dispersion(const Point& point, const Configurations&... configurations) const {
    const auto number_types(get_number_types());
    using return_t = decltype(matrix_times_vector_at_index(alpha_[0],
                                                           compute_dispersion(dispersion_[0], point, number_types, configurations...),
                                                           0));
    return_t return_value(0.);
    for(decltype(alpha_.size()) i(0); i < alpha_.size(); ++i) {
      if(!is_column_zero(alpha_[i], get_type(point))) {
        const auto short_range_dispersion(compute_dispersion(dispersion_[i], point, number_types, configurations...));
        return_value += matrix_times_vector_at_index(alpha_[i], short_range_dispersion, get_type(point));
      }
    }
    if(!is_column_zero(gamma_, get_type(point))) {
      const auto medium_range_dispersion(compute_dispersion(medium_range_dispersion_, point, number_types, configurations...));
      return_value += matrix_times_vector_at_index(gamma_, medium_range_dispersion, get_type(point));
    }
    return return_value;
  }

  template<typename Points, typename Configuration>
  auto compute_total_dispersion_vectorized(const Points& points, const Configuration& configuration, int nthreads) const {
    // Are any of the alpha/gamma rows non-zero?
    std::vector<bool> is_any_alpha_nonzero(alpha_.size());
    std::fill(is_any_alpha_nonzero.begin(), is_any_alpha_nonzero.end(), false);
    bool is_any_gamma_nonzero(false);
    for(typename Points::size_type i(0); i < size(points); ++i) {
      for(decltype(alpha_.size()) j(0); j < alpha_.size(); ++j) {
        if(!is_column_zero(alpha_[j], get_type(points[i]))) {
          is_any_alpha_nonzero[j] = true;
        }
      }
      if(!is_column_zero(gamma_, get_type(points[i]))) {
        is_any_gamma_nonzero = true;
      }
      if(std::all_of(is_any_alpha_nonzero.begin(),
                     is_any_alpha_nonzero.end(),
                     [](bool v) { return v; }) && is_any_gamma_nonzero) {
        break;
      }
    }

    // Precompute dispersion if required
    const auto number_types(get_number_types());
    using dispersion_t = decltype(compute_dispersion_for_fitting<true>(dispersion_[0], number_types, 1, configuration, points));
    std::vector<dispersion_t> short_range_dispersion{};
    dispersion_t medium_range_dispersion;
    for(decltype(dispersion_.size()) i(0); i < dispersion_.size(); ++i) {
      if(is_any_alpha_nonzero[i]) {
        short_range_dispersion.emplace_back(compute_dispersion_for_fitting<true>(dispersion_[i],
                                                                                 number_types,
                                                                                 nthreads,
                                                                                 configuration,
                                                                                 points));
      }
    }
    if(is_any_gamma_nonzero) {
      medium_range_dispersion = compute_dispersion_for_fitting<true>(medium_range_dispersion_, number_types, nthreads, configuration, points);
    }

    // Fill in the returned vector
    std::vector<decltype(matrix_times_vector_at_index(alpha_[0], short_range_dispersion[0][0], 0))> return_value(size(points));
    for(typename Points::size_type i(0); i < size(points); ++i) {
      // TODO: Put back commented code below??
      //if(!is_column_zero(alpha_, get_type(points[i])) || !is_column_zero(gamma_, get_type(points[i]))) {
        typename Configuration::size_type index_in_configuration(0);
        while(index_in_configuration < size(configuration) && !is_equal(configuration[index_in_configuration], points[i])) {
          ++index_in_configuration;
        }
        for(decltype(alpha_.size()) j(0); j < alpha_.size(); ++j) {
          if(!is_column_zero(alpha_[j], get_type(points[i]))) {
            return_value[i] += matrix_times_vector_at_index(alpha_[j],
                                                            short_range_dispersion[j][index_in_configuration < size(configuration) ? index_in_configuration : i + size(configuration)],
                                                            get_type(points[i]));
          }
        }
        if(!is_column_zero(gamma_, get_type(points[i]))) {
          return_value[i] += matrix_times_vector_at_index(gamma_, medium_range_dispersion[index_in_configuration < size(configuration) ? index_in_configuration : i + size(configuration)], get_type(points[i]));
        }
      //}
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
    using float_t = decltype(dispersion.get_maximum());
    if(dispersion.is_nonincreasing_after_lower_endpoint()) {
      return static_cast<float_t>(6.) * dispersion.get_maximum();
    } else {
      return std::numeric_limits<float_t>::infinity();
    }
  }
public:
  Truncated_exponential_family_model_over_window(const Window& window,
                                                 const Lambda& beta0,
                                                 Rcpp::List model,
                                                 Rcpp::CharacterVector medium_range_model,
                                                 Rcpp::List alpha,
                                                 Rcpp::NumericMatrix beta,
                                                 Rcpp::NumericMatrix gamma,
                                                 Rcpp::List covariates,
                                                 Rcpp::List short_range,
                                                 Rcpp::NumericMatrix medium_range,
                                                 Rcpp::NumericMatrix long_range,
                                                 unsigned long long int saturation):
    Model(beta0, model, medium_range_model, alpha, beta, gamma, covariates, short_range, medium_range, long_range, saturation),
    window_(window),
    beta_dot_covariates_maximum_(detail::compute_beta_dot_covariates_maximum(beta, Model::covariates_)),
    dot_dispersion_maximum_(positive_matrix_times_vector(Rcpp::as<Rcpp::NumericMatrix>(alpha[0]), get_dispersion_maximum(Model::dispersion_[0]))) {
    for(decltype(alpha.size()) i(1); i < alpha.size(); ++i) {
      const auto current_dot_dispersion_maximum(positive_matrix_times_vector(Rcpp::as<Rcpp::NumericMatrix>(alpha[i]), get_dispersion_maximum(Model::dispersion_[i])));
      std::transform(dot_dispersion_maximum_.begin(), dot_dispersion_maximum_.end(), current_dot_dispersion_maximum.begin(),
                     dot_dispersion_maximum_.begin(), std::plus<typename decltype(dot_dispersion_maximum_)::value_type>());
    }
    const auto gamma_dot_dispersion_maximum(positive_matrix_times_vector(gamma, get_dispersion_maximum(Model::medium_range_dispersion_)));
    std::transform(dot_dispersion_maximum_.begin(), dot_dispersion_maximum_.end(), gamma_dot_dispersion_maximum.begin(),
                   dot_dispersion_maximum_.begin(), std::plus<typename decltype(dot_dispersion_maximum_)::value_type>());
  }

  template<typename Point>
  auto get_log_normalized_bounding_intensity(const Point& point) const {
    auto beta_covariates_maximum(beta_dot_covariates_maximum_[get_type(point)]);
    auto beta_covariates(detail::compute_beta_dot_covariates(point, Model::beta_, Model::covariates_));
    return beta_covariates - beta_covariates_maximum;
  }

  // TODO: This should somehow be restricted to window.
  auto get_integral() const {
    double integral(0);
    const auto number_types(Model::get_number_types());
    using size_t = std::remove_cv_t<decltype(number_types)>;
    for(size_t type(0); type < number_types; ++type) {
      integral += Model::covariates_.get_integral_of_dot(window_, [beta0 = Model::beta0_[type], dot = dot_dispersion_maximum_[type]](auto x) {
        return std::exp(x + beta0 + dot);
      }, Model::beta_.get_row(type));
    }
    return integral;
  }

  template<typename Generator>
  auto sample_point_from_bounding_intensity(Generator& generator) const {
    // Compute exp(beta0)
    std::vector<double> exp_beta0(Model::beta0_.size());
    for(decltype(exp_beta0.size()) i(0); i < exp_beta0.size(); ++i) {
      exp_beta0[i] = std::exp(Model::beta0_[i]);
    }

    // Random distributions
    std::discrete_distribution<decltype(Model::beta0_.size())> random_type_distribution(exp_beta0.begin(), exp_beta0.end());
    std::exponential_distribution<double> exponential_distribution(1);

    // Sample type proportionally to the exp(beta0_i).
    const auto random_type(random_type_distribution(generator));

    while(true) {
      const auto sample(window_.sample(generator, random_type));
      if(exponential_distribution(generator) + get_log_normalized_bounding_intensity(sample) >= 0) {
        return sample;
      }
    }
  }

  auto get_upper_bound_approximate_ppp_intensity() const {
    const auto number_types(Model::get_number_types());
    std::vector<decltype(Model::beta0_[0] + beta_dot_covariates_maximum_[0])> upper_bound(number_types);
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
    std::vector<decltype(Model::beta0_[0] + beta_dot_covariates_maximum_[0])> upper_bound(number_types);
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
    if(std::all_of(Model::alpha_.begin(), Model::alpha_.end(), [&point](const auto alpha) {
          return(is_column_nonnegative(alpha, get_type(point)));
        }) && is_column_nonnegative(Model::gamma_, get_type(point))) {
      alpha_need_to_add = Model::compute_total_dispersion(point, l, l_complement) - dot_dispersion_maximum;
      if(alpha_need_to_add + exp_mark > 0.) {
        alpha_which_one_to_add_to = Model::compute_total_dispersion(point, l) - dot_dispersion_maximum;
        add_point(alpha_which_one_to_add_to + exp_mark > 0. ? l : l_complement, point);
      }
    } else if(std::all_of(Model::alpha_.begin(), Model::alpha_.end(), [&point](const auto alpha) {
          return(is_column_nonpositive(alpha, get_type(point)));
        }) && is_column_nonpositive(Model::gamma_, get_type(point))) {
      alpha_need_to_add = Model::compute_total_dispersion(point, l) - dot_dispersion_maximum;
      if(alpha_need_to_add + exp_mark > 0.) {
        alpha_which_one_to_add_to = Model::compute_total_dispersion(point, l, l_complement) - dot_dispersion_maximum;
        add_point(alpha_which_one_to_add_to + exp_mark > 0. ? l : l_complement, point);
      }
    } else {
      using dispersion_t = std::remove_cv_t<std::remove_reference_t<decltype(compute_dispersion(Model::dispersion_[0], point, Model::alpha_[0].nrow(), l))>>;
      std::vector<dispersion_t> short_range_dispersion_l, short_range_dispersion_u;
      dispersion_t medium_range_dispersion_l, medium_range_dispersion_u;
      alpha_need_to_add = -dot_dispersion_maximum;
      for(decltype(Model::alpha_.size()) i(0); i < Model::alpha_.size(); ++i) {
        if(!is_column_zero(Model::alpha_[i], get_type(point))) {
          // TODO: These two can definitely be computed faster by doing both at once.
          short_range_dispersion_u.emplace_back(compute_dispersion(Model::dispersion_[i], point, Model::alpha_[i].nrow(), l, l_complement));
          short_range_dispersion_l.emplace_back(compute_dispersion(Model::dispersion_[i], point, Model::alpha_[i].nrow(), l));
          alpha_need_to_add += positive_matrix_times_vector_at_index(Model::alpha_[i], short_range_dispersion_u[i], get_type(point));
          alpha_need_to_add += negative_matrix_times_vector_at_index(Model::alpha_[i], short_range_dispersion_l[i], get_type(point));
        }
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
        for(decltype(Model::alpha_.size()) i(0); i < Model::alpha_.size(); ++i) {
          if(!is_column_zero(Model::alpha_[i], get_type(point))) {
            alpha_which_one_to_add_to += positive_matrix_times_vector_at_index(Model::alpha_[i],
                                                                               short_range_dispersion_l[i], get_type(point));
            alpha_which_one_to_add_to += negative_matrix_times_vector_at_index(Model::alpha_[i],
                                                                               short_range_dispersion_u[i], get_type(point));
          }
        }
        if(!is_column_zero(Model::gamma_, get_type(point))) {
          alpha_which_one_to_add_to += positive_matrix_times_vector_at_index(Model::gamma_, medium_range_dispersion_l, get_type(point));
          alpha_which_one_to_add_to += negative_matrix_times_vector_at_index(Model::gamma_, medium_range_dispersion_u, get_type(point));
        }
        add_point(alpha_which_one_to_add_to + exp_mark > 0. ? l : l_complement, point);
      }
    }
    if(alpha_need_to_add > 0.) {
      Rcpp::stop(std::string("Internal error: Incorrect first alpha value, probably due to bad upper-bound to dispersion in the CFTP algorithm; alpha = ").append(std::to_string(alpha_need_to_add)).append(std::string(", point type = ")).append(std::to_string(get_type(point))).append(std::string(", upper bound = ")).append(std::to_string(dot_dispersion_maximum)));
    }
    if(alpha_which_one_to_add_to > 0.) {
      Rcpp::stop(std::string("Internal error: Incorrect second alpha value, probably due to bad upper-bound to dispersion in the CFTP algorithm; alpha = ").append(std::to_string(alpha_which_one_to_add_to)).append(std::string(", point type = ")).append(std::to_string(get_type(point))).append(std::string(", upper bound = ")).append(std::to_string(dot_dispersion_maximum)));
    }
  }

  auto compute_log_alpha_min_lower_bound(R_xlen_t type) const {
    using sum_t = decltype(std::abs(Model::alpha_[0](0, 0)));
    sum_t sum(0.);
    for(decltype(Model::alpha_.size()) i(0); i < Model::alpha_.size(); ++i) {
      sum_t sum_alpha(0.);
      for(decltype(Model::alpha_[i].ncol()) other_type(0); other_type < Model::alpha_[i].ncol(); ++other_type) {
        sum_alpha -= std::abs(Model::alpha_[i](type, other_type));
      }
      sum_alpha *= get_dispersion_maximum(Model::dispersion_[i]);
      sum += sum_alpha;
    }
    sum_t sum_gamma(0.);
    for(decltype(Model::gamma_.ncol()) other_type(0); other_type < Model::gamma_.ncol(); ++other_type) {
      sum_gamma -= std::abs(Model::gamma_(type, other_type));
    }
    return sum + sum_gamma * get_dispersion_maximum(Model::medium_range_dispersion_);
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

template<typename Configuration, typename Generator, typename Lambda>
struct approximate_draw_helper<Configuration, Generator, Truncated_exponential_family_model_over_window<Lambda>> {
  static auto get(Generator& generator, const Truncated_exponential_family_model_over_window<Lambda>& model) {
    const auto configuration(simulate_inhomogeneous_ppp<Configuration>(generator,
                                                                       model.get_window(),
                                                                       [&model](const auto& point) {
                                                                         return model.get_log_normalized_bounding_intensity(point);
                                                                       },
                                                                       model.get_upper_bound_approximate_ppp_intensity()));
    return configuration;
  }
};

} // namespace detail

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_EXPONENTIAL_FAMILY
