#ifndef INCLUDE_PPJSDM_PHI_DISPERSION
#define INCLUDE_PPJSDM_PHI_DISPERSION

#include <Rcpp.h>

#include <cmath> // std::exp, std::log
#include <limits> // std::numeric_limits
#include <type_traits> // std::remove_const, std::remove_cv
#include <utility> // std::forward
#include <vector> // std::vector

#include "compute_phi_distance.h"
#include "../configuration/configuration_manipulation.h"
#include "../configuration/get_number_points.h"
#include "../point/point_manipulation.h"

namespace ppjsdm {

template<typename Varphi>
class Varphi_model_papangelou: public Varphi {
public:
  template<typename... Args>
  Varphi_model_papangelou(Args&&... args): Varphi(std::forward<Args>(args)...) {}
  template<typename Configuration, typename Point>
  inline Rcpp::NumericVector compute(const Configuration& configuration, const Point& point, R_xlen_t number_types, R_xlen_t number_points) const {
    R_xlen_t same_type_points(0);
    for(R_xlen_t i(0); i < number_points; ++i) {
      const auto type_i(get_type(configuration[i]));
      if(type_i == get_type(point)) {
        ++same_type_points;
      }
    }

    // delta_dispersion and count_species are automatically 0-initialized
    Rcpp::NumericVector delta_dispersion(number_types);
    Rcpp::IntegerVector count_species(number_types);

    for(R_xlen_t i(0); i < number_points; ++i) {
      const auto point_i(configuration[i]);
      const auto type_i(get_type(point_i));
      ++count_species[type_i];

      double average_dispersion(0);

      for(R_xlen_t j(0); j < number_points; ++j) {
        const auto point_j(configuration[j]);
        const auto type_j(get_type(point_j));
        if(type_j == get_type(point) && j != i) {
          const auto phi_distance(compute_phi_distance(point_i, point_j, *this));

          average_dispersion += phi_distance;
        }
      }
      if(type_i == get_type(point)) {
        if(same_type_points > 1) {
          average_dispersion /= same_type_points - 1;
        }
      } else {
        if(same_type_points > 0) {
          average_dispersion /= same_type_points;
        }
      }
      delta_dispersion[type_i] += average_dispersion - compute_phi_distance(point_i, point, *this);
    }

    for(R_xlen_t i(0); i < number_types; ++i) {
      const auto count_species_i(count_species[i]);
      if(count_species_i > 0) {
        if(i == get_type(point)) {
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
private:
  auto access_square_radii(R_xlen_t i, R_xlen_t j) const {
    return square_radii_[i * dim_ + j];
  }
  void set_square_radii(R_xlen_t i, R_xlen_t j, double r) {
    square_radii_[i * dim_ + j] = r * r;
  }
public:
  Geyer_papangelou(Rcpp::NumericMatrix radius, double saturation): dim_(radius.ncol()), square_radii_(dim_ * dim_), saturation_(saturation) {
    for(R_xlen_t i(0); i < dim_; ++i) {
      for(R_xlen_t j(0); j < dim_; ++j) {
        set_square_radii(i, j, radius(i, j));
      }
    }
  }
  template<typename Configuration, typename Point>
  inline Rcpp::NumericVector compute(const Configuration& configuration, const Point& point, R_xlen_t number_types, R_xlen_t number_points) const {
    Rcpp::NumericVector delta_dispersion(number_types);
    const auto type_point(get_type(point));
    for(R_xlen_t i(0); i < number_points; ++i) {
      const auto point_i(configuration[i]);
      const auto type_i(get_type(point_i));
      if(delta_dispersion[type_i] < saturation_) {
        const auto delta_x(get_x(point_i) - get_x(point));
        const auto delta_y(get_y(point_i) - get_y(point));
        const auto square_distance(delta_x * delta_x + delta_y * delta_y);
        if(square_distance <= access_square_radii(type_i, type_point)) {
          delta_dispersion[type_i] += 1.;
        }
      }
    }
    return delta_dispersion;
  }

private:
  R_xlen_t dim_;
  std::vector<double> square_radii_;
  double saturation_;
};

template<typename Varphi>
class Nearest_neighbour_papangelou: public Varphi {
public:
  template<typename... Args>
  Nearest_neighbour_papangelou(Args&&... args): Varphi(std::forward<Args>(args)...) {}

  template<typename Configuration, typename Point>
  inline Rcpp::NumericVector compute(const Configuration& configuration, const Point& point, R_xlen_t number_types, R_xlen_t number_points) const {
    Rcpp::NumericVector delta_dispersion(Rcpp::no_init(number_types));

    auto configuration_plus(configuration);
    add_point(configuration_plus, point);

    for(R_xlen_t i(0); i < number_types; ++i) {
      if(i == get_type(point)) {
        const double d(-compute_dispersion(configuration, i, i, number_points));
        delta_dispersion[i] = d + compute_dispersion(configuration_plus, i, i, number_points + 1);
      } else {
        auto d(-compute_dispersion(configuration, i, get_type(point), number_points));
        d -= compute_dispersion(configuration, get_type(point), i, number_points);
        delta_dispersion[i] = 0.5 * (d + compute_dispersion(configuration_plus, i, get_type(point), number_points + 1) + compute_dispersion(configuration_plus, get_type(point), i, number_points + 1));
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
      const auto point_i(configuration[i]);
      const auto type_i(get_type(point_i));
      if(type_i == k1) {
        ++count_species;
        double min_distance(std::numeric_limits<double>::max());
        for(R_xlen_t j(0); j < number_points; ++j) {
          if(j != i) {
            const auto point_j(configuration[j]);
            const auto type_j(get_type(point_j));
            if(type_j == k2) {
              const auto d(compute_phi_distance(point_i, point_j, *this));
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

template<typename Varphi, typename Lambda, typename Nu, typename Alpha>
class Exponential_family_model: public Varphi {
public:
  template<typename... Args>
  Exponential_family_model(const Lambda& lambda, const Nu& nu, const Alpha& alpha, Args&&... args): Varphi(std::forward<Args>(args)...), lambda_(lambda), nu_(nu), alpha_(alpha) {}
  Exponential_family_model(const Lambda& lambda, const Nu& nu, const Alpha& alpha, const Varphi& varphi): Varphi(varphi), lambda_(lambda), nu_(nu), alpha_(alpha) {}
  Exponential_family_model(const Lambda& lambda, const Nu& nu, const Alpha& alpha, Varphi&& varphi): Varphi(std::move(varphi)), lambda_(lambda), nu_(nu), alpha_(alpha) {}

  template<typename Configuration, typename Point>
  double compute_papangelou(const Configuration& configuration, const Point& point, R_xlen_t number_types) const {
    const auto delta_D(Varphi::compute(configuration, point, number_types, size(configuration)));

    double inner_product(0);
    const auto type(get_type(point));
    for(R_xlen_t i(0); i < number_types; ++i) {
      inner_product += alpha_(i, type) * delta_D[i];
    }
    return lambda_[type] * std::exp(inner_product + (1 - nu_[type]) * std::log(get_number_points(configuration, type) + 1));
  }
private:
  Lambda lambda_;
  Nu nu_;
  Alpha alpha_;
};

const constexpr char* const models[] = {
  "identity",
  "Strauss",
  "Geyer",
  "neighbour"
};

template<typename F>
inline auto call_on_papangelou(Rcpp::CharacterVector model, Rcpp::NumericMatrix radius, const F& f) {
  const auto model_string(model[0]);
  if(model_string == models[0]) {
    return f(Varphi_model_papangelou<varphi::Identity>{});
  } else if(model_string == models[1]) {
    return f(Varphi_model_papangelou<varphi::Strauss>(radius));
  } else if(model_string == models[2]) {
    return f(Geyer_papangelou(radius, 2.0));
  } else if(model_string == models[3]) {
    return f(Nearest_neighbour_papangelou<varphi::Identity>{});
  } else {
    Rcpp::stop("Incorrect model entered.\n");
  }
}

template<typename F, typename L, typename N, typename... Args>
inline auto call_on_model(Rcpp::CharacterVector model, Rcpp::NumericMatrix alpha, const L& lambda, const N& nu, Rcpp::NumericMatrix radius, const F& f) {
  return call_on_papangelou(model, radius, [&alpha, &lambda, &nu, &f](auto&& varphi){
    const Exponential_family_model<std::remove_cv_t<std::remove_reference_t<decltype(varphi)>>, L, N, Rcpp::NumericMatrix> model(lambda, nu, alpha, std::forward<decltype(varphi)>(varphi));
    return f(model);
  });
}


} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_PHI_DISPERSION
