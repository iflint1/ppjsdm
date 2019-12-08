#ifndef INCLUDE_PPJSDM_PHI_DISPERSION
#define INCLUDE_PPJSDM_PHI_DISPERSION

#include <Rcpp.h>

#include "compute_phi_distance.h"
// #include "configuration_utilities.h"

template<typename V, typename X, typename Y, typename T>
[[nodiscard]] inline Rcpp::NumericMatrix compute_phi_dispersion_helper(const X& x, const Y& y, const T& types_vector, int number_types, const V& varphi) {
  // TODO: Might be able to avoid recomputing this every time, marginal efficiency gain.
  const R_xlen_t number_points{x.size()};

  // dispersion and count_species are automatically 0-initialized
  Rcpp::NumericMatrix dispersion(number_types, number_types);
  Rcpp::IntegerVector count_species(number_types);

  for(R_xlen_t i{0}; i < number_points; ++i) {
    const auto type_i{types_vector[i] - 1};
    ++count_species[type_i];

    for(R_xlen_t j{i + 1}; j < number_points; ++j) {
      const auto type_j{types_vector[j] - 1};

      const auto phi_distance{compute_phi_distance(x[i], y[i], x[j], y[j], varphi)};

      dispersion(type_i, type_j) += phi_distance;
      if(type_i != type_j) {
        dispersion(type_j, type_i) += phi_distance;
      }
    }
  }

  for(int i{0}; i < number_types; ++i) {
    const auto count_species_i{count_species[i]};
    if(count_species_i > 0) {
      if(count_species_i != 1) {
        dispersion(i, i) *= 2. / (count_species_i * (count_species_i - 1));
      }
      for(int j{i + 1}; j < number_types; ++j) {
        const auto count_species_j{count_species[j]};
        if(count_species_j > 0) {
          dispersion(i, j) /= count_species_i * count_species_j;
          dispersion(j, i) /= count_species_i * count_species_j;
        }
      }
    }
  }

  return dispersion;
}

template<typename V>
[[nodiscard]] inline Rcpp::NumericMatrix compute_phi_dispersion_helper(Rcpp::S4 configuration, int number_types, const V& varphi) {
  return compute_phi_dispersion_helper(Rcpp::NumericVector(configuration.slot("x")),
                                       Rcpp::NumericVector(configuration.slot("y")),
                                       Rcpp::IntegerVector(configuration.slot("types")),
                                       number_types,
                                       varphi);
}

template<typename V, typename X, typename Y, typename T>
[[nodiscard]] inline Rcpp::NumericVector compute_delta_phi_dispersion_helper(const X& x, const Y& y, const T& types_vector, Rcpp::NumericVector location, R_xlen_t type, R_xlen_t number_types, const V& varphi) {
  // TODO: Might be able to avoid recomputing this every time, marginal efficiency gain.
  const auto number_points{x.size()};

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
        const auto phi_distance{compute_phi_distance(x[i], y[i], x[j], y[j], varphi)};

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
    delta_dispersion[type_i] += average_dispersion - compute_phi_distance(x[i], y[i], location[0], location[1], varphi);
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

template<typename V>
[[nodiscard]] inline Rcpp::NumericVector compute_delta_phi_dispersion_helper(Rcpp::S4 configuration, Rcpp::NumericVector location, R_xlen_t type, R_xlen_t number_types, const V& varphi) {
  return compute_delta_phi_dispersion_helper(Rcpp::NumericVector(configuration.slot("x")),
                                          Rcpp::NumericVector(configuration.slot("y")),
                                          Rcpp::IntegerVector(configuration.slot("types")),
                                          location,
                                          type,
                                          number_types,
                                          varphi);
}

#endif // INCLUDE_PPJSDM_PHI_DISPERSION
