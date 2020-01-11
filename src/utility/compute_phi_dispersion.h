#ifndef INCLUDE_PPJSDM_PHI_DISPERSION
#define INCLUDE_PPJSDM_PHI_DISPERSION

#include <Rcpp.h>

#include <cmath> // std::exp, std::log
#include <limits> // std::numeric_limits
#include <type_traits> // std::remove_const, std::remove_cv
#include <utility> // std::forward

#include "compute_phi_distance.h"

namespace ppjsdm {

template<typename Varphi>
class Varphi_model_papangelou: public Varphi {
public:
  template<typename... Args>
  Varphi_model_papangelou(Args&&... args): Varphi(std::forward<Args>(args)...) {}
  template<typename Configuration>
  inline Rcpp::NumericVector compute(const Configuration& configuration, Rcpp::NumericVector location, R_xlen_t type, R_xlen_t number_types, R_xlen_t number_points) const {
    R_xlen_t same_type_points(0);
    for(R_xlen_t i(0); i < number_points; ++i) {
      const auto type_i(configuration.types(i) - 1);
      if(type_i == type) {
        ++same_type_points;
      }
    }

    // delta_dispersion and count_species are automatically 0-initialized
    Rcpp::NumericVector delta_dispersion(number_types);
    Rcpp::IntegerVector count_species(number_types);

    for(R_xlen_t i(0); i < number_points; ++i) {
      const auto type_i(configuration.types(i) - 1);
      ++count_species[type_i];

      double average_dispersion(0);

      for(R_xlen_t j(0); j < number_points; ++j) {
        const auto type_j(configuration.types(j) - 1);
        if(type_j == type && j != i) {
          const auto phi_distance(compute_phi_distance(configuration.x(i), configuration.y(i), configuration.x(j), configuration.y(j), *this));

          average_dispersion += phi_distance;
        }
      }
      if(type_i == type) {
        if(same_type_points > 1) {
          average_dispersion /= same_type_points - 1;
        }
      } else {
        if(same_type_points > 0) {
          average_dispersion /= same_type_points;
        }
      }
      delta_dispersion[type_i] += average_dispersion - compute_phi_distance(configuration.x(i), configuration.y(i), location[0], location[1], *this);
    }

    for(R_xlen_t i(0); i < number_types; ++i) {
      const auto count_species_i(count_species[i]);
      if(count_species_i > 0) {
        if(i == type) {
          delta_dispersion[i] *= 2. / (count_species_i * (count_species_i + 1));
        } else {
          delta_dispersion[i] /= (count_species_i * (same_type_points + 1));
        }
      }
    }

    return delta_dispersion;
  }
};

class Geyer_papangelou {
public:
  Geyer_papangelou(double radius, double saturation) noexcept: square_radius_(radius * radius), saturation_(saturation) {}
  template<typename Configuration>
  inline Rcpp::NumericVector compute(const Configuration& configuration, Rcpp::NumericVector location,R_xlen_t, R_xlen_t number_types, R_xlen_t number_points) const {
    Rcpp::NumericVector delta_dispersion(number_types);
    for(R_xlen_t i(0); i < number_points; ++i) {
      const auto type_i(configuration.types(i) - 1);
      if(delta_dispersion[type_i] < saturation_) {
        const auto delta_x(configuration.x(i) - location[0]);
        const auto delta_y(configuration.y(i) - location[1]);
        const auto square_distance(delta_x * delta_x + delta_y * delta_y);
        if(square_distance <= square_radius_) {
          delta_dispersion[type_i] += 1.;
        }
      }
    }

    return delta_dispersion;
  }

private:
  double square_radius_;
  double saturation_;
};

template<typename Varphi>
class Nearest_neighbour_papangelou: public Varphi {
public:
  template<typename... Args>
  Nearest_neighbour_papangelou(Args&&... args): Varphi(std::forward<Args>(args)...) {}

  template<typename Configuration>
  inline Rcpp::NumericVector compute(const Configuration& configuration, Rcpp::NumericVector location, R_xlen_t type, R_xlen_t number_types, R_xlen_t number_points) const {
    Rcpp::NumericVector delta_dispersion(Rcpp::no_init(number_types));

    auto configuration_plus(configuration);
    emplace_point(configuration_plus, location[0], location[1], type + 1);

    for(R_xlen_t i(0); i < number_types; ++i) {
      if(i == type) {
        const double d(-compute_dispersion(configuration, i, i, number_points));
        delta_dispersion[i] = d + compute_dispersion(configuration_plus, i, i, number_points + 1);
      } else {
        auto d(-compute_dispersion(configuration, i, type, number_points));
        d -= compute_dispersion(configuration, type, i, number_points);
        delta_dispersion[i] = 0.5 * (d + compute_dispersion(configuration_plus, i, type, number_points + 1) + compute_dispersion(configuration_plus, type, i, number_points + 1));
      }
    }
    return delta_dispersion;
  }
private:
  template<typename Configuration>
  inline double compute_dispersion(const Configuration& configuration,
                                                 R_xlen_t k1, R_xlen_t k2,
                                                 R_xlen_t number_points) const {
    int count_species(0);
    double dispersion(0);

    for(R_xlen_t i(0); i < number_points; ++i) {
      const auto type_i(configuration.types(i) - 1);
      if(type_i == k1) {
        ++count_species;
        double min_distance(std::numeric_limits<double>::max());
        for(R_xlen_t j(0); j < number_points; ++j) {
          if(j != i) {
            const auto type_j(configuration.types(j) - 1);
            if(type_j == k2) {
              const auto d(compute_phi_distance(configuration.x(i), configuration.y(i), configuration.x(j), configuration.y(j), *this));
              if(min_distance > d) {
                min_distance = d;
              }
            }
          }
        }
        if(min_distance < std::numeric_limits<double>::max()) {
          dispersion += min_distance;
        }
      }
    }

    if(count_species > 0) {
      dispersion /= static_cast<double>(count_species);
    }
    return dispersion;
  }
};

template<typename Varphi, typename U, typename V, typename W>
class Exponential_family_model: public Varphi {
public:
  template<typename... Args>
  Exponential_family_model(const U& lambda, const V& nu, const W& alpha, Args&&... args): Varphi(std::forward<Args>(args)...), lambda_(lambda), nu_(nu), alpha_(alpha) {}
  Exponential_family_model(const U& lambda, const V& nu, const W& alpha, const Varphi& varphi): Varphi(varphi), lambda_(lambda), nu_(nu), alpha_(alpha) {}
  Exponential_family_model(const U& lambda, const V& nu, const W& alpha, Varphi&& varphi): Varphi(std::move(varphi)), lambda_(lambda), nu_(nu), alpha_(alpha) {}

  template<typename Configuration>
  double compute_papangelou(const Configuration& configuration, Rcpp::NumericVector location, R_xlen_t type, R_xlen_t number_types) const {
    const auto delta_D(Varphi::compute(configuration, location, type, number_types, size(configuration)));

    double inner_product(0);
    for(R_xlen_t i(0); i < number_types; ++i) {
      inner_product += alpha_(i, type) * delta_D[i];
    }
    // TODO: replace x.size() with "number of points of type `type`".
    // TODO: Faster with n ^ nu?
    return lambda_[type] * std::exp(inner_product + (1 - nu_[type]) * std::log(size(configuration) + 1));
  }
private:
  U lambda_;
  V nu_;
  W alpha_;
};

template<typename F>
inline auto call_on_papangelou(Rcpp::CharacterVector model, double radius, const F& f) {
  const auto model_string(model[0]);
  if(model_string == "identity") {
    return f(Varphi_model_papangelou<varphi::Identity>{});
  } else if(model_string == "Strauss") {
    return f(Varphi_model_papangelou<varphi::Strauss>(radius));
  } else if(model_string == "Geyer") {
    return f(Geyer_papangelou(radius, 2.0));
  } else if(model_string == "neighbour") {
    return f(Nearest_neighbour_papangelou<varphi::Identity>{});
  } else {
    Rcpp::stop("Incorrect model entered.\n");
  }
}

template<typename F, typename L, typename N, typename... Args>
inline auto call_on_model(Rcpp::CharacterVector model, Rcpp::NumericMatrix alpha, const L& lambda, const N& nu, double radius, const F& f) {
  return call_on_papangelou(model, radius, [&alpha, &lambda, &nu, &f](auto&& varphi){
    const Exponential_family_model<std::remove_cv_t<std::remove_reference_t<decltype(varphi)>>, L, N, Rcpp::NumericMatrix> model(lambda, nu, alpha, std::forward<decltype(varphi)>(varphi));
    return f(model);
  });
}


} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_PHI_DISPERSION
