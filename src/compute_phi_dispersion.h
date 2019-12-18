#ifndef INCLUDE_PPJSDM_PHI_DISPERSION
#define INCLUDE_PPJSDM_PHI_DISPERSION

#include <Rcpp.h>

#include <utility> // std::forward

#include "compute_phi_distance.h"
// #include "configuration_utilities.h"

// template<typename V, typename X, typename Y, typename T>
// [[nodiscard]] inline Rcpp::NumericMatrix compute_phi_dispersion_helper(const X& x, const Y& y, const T& types_vector, int number_types, const V& varphi) {
//   // TODO: Might be able to avoid recomputing this every time, marginal efficiency gain.
//   const R_xlen_t number_points{x.size()};
//
//   // dispersion and count_species are automatically 0-initialized
//   Rcpp::NumericMatrix dispersion(number_types, number_types);
//   Rcpp::IntegerVector count_species(number_types);
//
//   for(R_xlen_t i{0}; i < number_points; ++i) {
//     const auto type_i{types_vector[i] - 1};
//     ++count_species[type_i];
//
//     for(R_xlen_t j{i + 1}; j < number_points; ++j) {
//       const auto type_j{types_vector[j] - 1};
//
//       const auto phi_distance{compute_phi_distance(x[i], y[i], x[j], y[j], varphi)};
//
//       dispersion(type_i, type_j) += phi_distance;
//       if(type_i != type_j) {
//         dispersion(type_j, type_i) += phi_distance;
//       }
//     }
//   }
//
//   for(int i{0}; i < number_types; ++i) {
//     const auto count_species_i{count_species[i]};
//     if(count_species_i > 0) {
//       if(count_species_i != 1) {
//         dispersion(i, i) *= 2. / (count_species_i * (count_species_i - 1));
//       }
//       for(int j{i + 1}; j < number_types; ++j) {
//         const auto count_species_j{count_species[j]};
//         if(count_species_j > 0) {
//           dispersion(i, j) /= count_species_i * count_species_j;
//           dispersion(j, i) /= count_species_i * count_species_j;
//         }
//       }
//     }
//   }
//
//   return dispersion;
// }
//
// template<typename V>
// [[nodiscard]] inline Rcpp::NumericMatrix compute_phi_dispersion_helper(Rcpp::S4 configuration, int number_types, const V& varphi) {
//   return compute_phi_dispersion_helper(Rcpp::NumericVector(configuration.slot("x")),
//                                        Rcpp::NumericVector(configuration.slot("y")),
//                                        Rcpp::IntegerVector(configuration.slot("types")),
//                                        number_types,
//                                        varphi);
// }

template<typename Varphi>
class Varphi_model_papangelou: public Varphi {
public:
  template<typename... Args>
  Varphi_model_papangelou(Args&&... args): Varphi{std::forward<Args>(args)...} {}
  template<typename X, typename Y, typename T>
  [[nodiscard]] inline Rcpp::NumericVector compute(const X& x, const Y& y, const T& types_vector, Rcpp::NumericVector location, R_xlen_t type, R_xlen_t number_types) const {
    // TODO: Might be able to avoid recomputing this every time, marginal efficiency gain.
    const auto number_points = R_xlen_t(x.size());

    R_xlen_t same_type_points{0};
    for(R_xlen_t i{0}; i < number_points; ++i) {
      const auto type_i{types_vector[i] - 1};
      if(type_i == type) {
        ++same_type_points;
      }
    }

    // delta_dispersion and count_species are automatically 0-initialized
    Rcpp::NumericVector delta_dispersion(number_types);
    Rcpp::IntegerVector count_species(number_types);

    for(R_xlen_t i{0}; i < number_points; ++i) {
      const auto type_i{types_vector[i] - 1};
      ++count_species[type_i];

      double average_dispersion{0};

      for(R_xlen_t j{0}; j < number_points; ++j) {
        const auto type_j{types_vector[j] - 1};
        if(type_j == type && j != i) {
          const auto phi_distance{compute_phi_distance(x[i], y[i], x[j], y[j], *this)};

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
      delta_dispersion[type_i] += average_dispersion - compute_phi_distance(x[i], y[i], location[0], location[1], *this);
    }

    for(R_xlen_t i{0}; i < number_types; ++i) {
      const auto count_species_i{count_species[i]};
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

// template<typename V>
// [[nodiscard]] inline Rcpp::NumericVector compute_delta_phi_dispersion_helper(Rcpp::S4 configuration, Rcpp::NumericVector location, R_xlen_t type, R_xlen_t number_types, const V& varphi) {
//   return compute_delta_phi_dispersion_helper(Rcpp::NumericVector(configuration.slot("x")),
//                                           Rcpp::NumericVector(configuration.slot("y")),
//                                           Rcpp::IntegerVector(configuration.slot("types")),
//                                           location,
//                                           type,
//                                           number_types,
//                                           varphi);
// }

class Geyer_papangelou {
public:
  Geyer_papangelou(double radius, double saturation) noexcept: square_radius_{radius * radius}, saturation_{saturation} {}
  template<typename X, typename Y, typename T>
  [[nodiscard]] inline Rcpp::NumericVector compute(const X& x, const Y& y, const T& types_vector, Rcpp::NumericVector location, R_xlen_t type, R_xlen_t number_types) const {
    // TODO: Might be able to avoid recomputing this every time, marginal efficiency gain.
    const auto number_points = R_xlen_t(x.size());

    Rcpp::NumericVector delta_dispersion(number_types);
    for(R_xlen_t i{0}; i < number_points; ++i) {
      const auto type_i{types_vector[i] - 1};
      if(delta_dispersion[type_i] < saturation_) {
        const auto delta_x{x[i] - location[0]};
        const auto delta_y{y[i] - location[1]};
        const auto square_distance{delta_x * delta_x + delta_y * delta_y};
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
  Nearest_neighbour_papangelou(Args&&... args): Varphi{std::forward<Args>(args)...} {}

  template<typename X, typename Y, typename T>
  [[nodiscard]] inline Rcpp::NumericVector compute(const X& x, const Y& y, const T& types_vector, Rcpp::NumericVector location, R_xlen_t type, R_xlen_t number_types) const {
    // TODO: Might want to not initialize
    Rcpp::NumericVector delta_dispersion(number_types);

    for(R_xlen_t i{0}; i < number_types; ++i) {
      if(i == type) {
        const auto d{-compute_dispersion(x, y, types_vector, i, i, number_types)};
        auto copy_x{x};
        auto copy_y{y};
        auto copy_types{types_vector};
        copy_x.push_back(location[0]);
        copy_y.push_back(location[1]);
        copy_types.push_back(type + 1);

        delta_dispersion[i] = d + compute_dispersion(copy_x, copy_y, copy_types, i, i, number_types);
      } else {
        auto d{-compute_dispersion(x, y, types_vector, i, type, number_types)};
        d -= compute_dispersion(x, y, types_vector, type, i, number_types);
        auto copy_x{x};
        auto copy_y{y};
        auto copy_types{types_vector};
        copy_x.push_back(location[0]);
        copy_y.push_back(location[1]);
        copy_types.push_back(type + 1);

        delta_dispersion[i] = 0.5 * (d + compute_dispersion(copy_x, copy_y, copy_types, i, type, number_types) + compute_dispersion(copy_x, copy_y, copy_types, type, i, number_types));
      }
    }

    return delta_dispersion;
  }
private:
  template<typename X, typename Y, typename T>
  [[nodiscard]] inline double compute_dispersion(const X& x, const Y& y,
                                                 const T& types_vector,
                                                 R_xlen_t k1, R_xlen_t k2,
                                                 R_xlen_t number_types) const {
    // TODO: Might be able to avoid recomputing this every time, marginal efficiency gain.
    const auto number_points = R_xlen_t(x.size());

    // count_species automatically 0-initialized
    R_xlen_t count_species;
    double dispersion{0};

    for(R_xlen_t i{0}; i < number_points; ++i) {
      const auto type_i{types_vector[i] - 1};
      if(type_i == k1) {
        ++count_species;
        // TODO: Change
        double min_distance{999999};
        for(R_xlen_t j{0}; j < number_points; ++j) {
          if(j != i) {
            const auto type_j{types_vector[j] - 1};
            if(type_j == k2) {
              const auto d{compute_phi_distance(x[i], y[i], x[j], y[j], *this)};
              if(min_distance > d) {
                min_distance = d;
              }
            }
          }
        }
        dispersion += min_distance;
      }
    }

    dispersion /= static_cast<double>(count_species);
  }
};

// [[nodiscard]] inline Rcpp::NumericVector compute_delta_phi_dispersion_geyer_helper(Rcpp::S4 configuration, Rcpp::NumericVector location, R_xlen_t type, R_xlen_t number_types, double square_radius, double saturation) {
//   return compute_delta_phi_dispersion_geyer_helper(Rcpp::NumericVector(configuration.slot("x")),
//                                              Rcpp::NumericVector(configuration.slot("y")),
//                                              Rcpp::IntegerVector(configuration.slot("types")),
//                                              location,
//                                              type,
//                                              number_types,
//                                              square_radius,
//                                              saturation);
// }

template<typename Varphi, typename U, typename V>
class Exponential_family_model: public Varphi {
public:
  template<typename... Args>
  Exponential_family_model(const U& lambda, const V& alpha, Args&&... args): Varphi(std::forward<Args>(args)...), lambda_{lambda}, alpha_{alpha} {}
  template<typename X, typename Y, typename T>
  [[nodiscard]] double compute_papangelou(const X& x, const Y& y, const T& types_vector, Rcpp::NumericVector location, R_xlen_t type, R_xlen_t number_types) const {
    const auto delta_D{Varphi::compute(x, y, types_vector, location, type, number_types)};

    double inner_product{0};
    for(R_xlen_t i{0}; i < number_types; ++i) {
      inner_product += alpha_(i, type) * delta_D[i];
    }

    return lambda_[type] * std::exp(inner_product);
  }
private:
  U lambda_;
  V alpha_;
};

#endif // INCLUDE_PPJSDM_PHI_DISPERSION
