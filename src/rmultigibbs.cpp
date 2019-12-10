#include <Rcpp.h>
#include <Rmath.h>

#include <cmath> // std::exp
#include <vector> // std::vector
#include <cstring> // std::strcmp

#include "compute_phi_dispersion.h"
#include "configuration_utilities.h"
#include "rbinomialpp_single.h"
#include "make_default_types.h"
#include "window_utilities.h"

// template<typename X, typename Y, typename T, typename V>
// [[nodiscard]] double compute_papangelou(const X& x, const Y& y, const T& types_vector, Rcpp::NumericVector location, R_xlen_t type, Rcpp::NumericMatrix alpha, Rcpp::NumericVector lambda, R_xlen_t number_types, const V& varphi) {
//
//   const auto delta_D{compute_delta_phi_dispersion_helper(x, y, types_vector, location, type, number_types, varphi)};
//
//   double inner_product{0};
//   for(R_xlen_t i{0}; i < number_types; ++i) {
//     inner_product += alpha(i, type) * delta_D[i];
//   }
//
//   return lambda[type] * std::exp(inner_product);
// }

// template<typename X, typename Y, typename T>
// [[nodiscard]] double compute_papangelou_geyer(const X& x, const Y& y, const T& types_vector, Rcpp::NumericVector location, R_xlen_t type, Rcpp::NumericMatrix alpha, Rcpp::NumericVector lambda, double saturation, R_xlen_t number_types, double square_radius) {
//
//   const auto delta_D{compute_delta_phi_dispersion_geyer_helper(x, y, types_vector, location, type, number_types, square_radius, saturation)};
//
//   double inner_product{0};
//   for(R_xlen_t i{0}; i < number_types; ++i) {
//     inner_product += alpha(i, type) *  delta_D[i];
//   }
//
//   return lambda[type] * std::exp(inner_product);
// }

template<typename V, typename S>
[[nodiscard]] inline SEXP rmultigibbs_helper2(const S& window, Rcpp::NumericMatrix alpha, Rcpp::NumericVector lambda, R_xlen_t steps, R_xlen_t nsim, Rcpp::Nullable<Rcpp::CharacterVector> types, bool drop, const V& varphi) {
  const auto number_types{lambda.size()};
  const auto volume{window.volume()};
  const auto prob{0.5};

  if(types.isNull()) {
    types = Rcpp::wrap(make_default_types(number_types));
  }
  const auto types_char_vector{types.as()};

  Rcpp::List samples{nsim};
  // TODO: This can likely be made faster by using one unique vector, since all 3 have the same length.
  // TODO: Make an std_vector_configuration wrapper to make everything more straightforward.
  std::vector<double> x;
  std::vector<double> y;
  std::vector<int> types_vector;

  for(R_xlen_t i{0}; i < nsim; ++i) {
    R_xlen_t total_number{0};

    for(R_xlen_t step{0}; step < steps; ++step) {
      const auto u{unif_rand()};
      const auto v{unif_rand()};
      if(u <= prob) {
        const R_xlen_t type{Rcpp::sample(number_types, 1, false, R_NilValue, false)[0]};
        const auto location_pair{window.sample()};
        const Rcpp::NumericVector location{location_pair.first, location_pair.second};

        // TODO: You can avoid taking the exp by reorganising the ratio, and sampling an exponential r.v. instead.
        const auto papangelou{varphi.compute_papangelou(x, y, types_vector, location, type, number_types)};
        //const auto papangelou{compute_papangelou_geyer(x, y, types_vector, location, type, alpha, lambda, 2., number_types, 0.02 * 0.02)};
        const auto birth_ratio{papangelou * (1 - prob) * volume * number_types / (prob * (1 + total_number))};

        if(v <= birth_ratio) {
          // Add point to configuration
          x.push_back(location[0]);
          y.push_back(location[1]);
          types_vector.push_back(type + 1);
          ++total_number;
        }
      } else {
        if(total_number != 0) {
          const R_xlen_t index{Rcpp::sample(x.size(), 1, false, R_NilValue, false)[0]};
          const Rcpp::NumericVector saved_location{x[index], y[index]};
          const R_xlen_t saved_type{types_vector[index] - 1};
          x.erase(x.begin() + index);
          y.erase(y.begin() + index);
          types_vector.erase(types_vector.begin() + index);

          const auto papangelou{varphi.compute_papangelou(x, y, types_vector, saved_location, saved_type, number_types)};
          //const auto papangelou{compute_papangelou_geyer(x, y, types_vector, saved_location, saved_type, alpha, lambda, 2., number_types, 0.02 * 0.02)};
          const auto death_ratio{prob * total_number / (number_types * (1 - prob) * volume * papangelou)};
          if(v <= death_ratio) {
            --total_number;
          } else {
            x.push_back(saved_location[0]);
            y.push_back(saved_location[1]);
            types_vector.push_back(saved_type + 1);
          }
        }
      }
    }

    samples[i] = make_configuration(Rcpp::wrap(x), Rcpp::wrap(y), Rcpp::wrap(types_vector), types_char_vector);
  }

  if(nsim == 1 && drop == true) {
    return Rcpp::wrap(samples[0]);
  } else {
    return Rcpp::wrap(samples);
  }
}

template<Window WindowType>
[[nodiscard]] inline SEXP rmultigibbs_helper(SEXP window, Rcpp::NumericMatrix alpha, Rcpp::NumericVector lambda, double radius, R_xlen_t steps, R_xlen_t nsim, Rcpp::Nullable<Rcpp::CharacterVector> types, Rcpp::CharacterVector model, bool drop) {
  const Window_wrapper<WindowType> window_wrapper{window};
  // TODO: Better string comparison below
  if(model[0] == "identity") {
    // TODO: Might use C++17 class type deduction
    const auto varphi{Exponential_family_model<Varphi_model_papangelou<varphi::Identity>, decltype(lambda), decltype(alpha)>{lambda, alpha}};
    return rmultigibbs_helper2(window_wrapper, alpha, lambda, steps, nsim, types, drop, varphi);
  } else if(model[0] == "Strauss") {
    const auto varphi{Exponential_family_model<Varphi_model_papangelou<varphi::Strauss>, decltype(lambda), decltype(alpha)>{lambda, alpha, radius * radius}};
    return rmultigibbs_helper2(window_wrapper, alpha, lambda, steps, nsim, types, drop, varphi);
  } else if(model[0] == "Geyer") {
    const auto varphi{Exponential_family_model<Geyer_papangelou, decltype(lambda), decltype(alpha)>{lambda, alpha, radius * radius, 2.0}};
    return rmultigibbs_helper2(window_wrapper, alpha, lambda, steps, nsim, types, drop, varphi);
  } else {
    Rcpp::stop("Incorrect model entered.\n");
  }
}

//' Sample a multivariate Gibbs point processes
//'
//' @param window The window.
//' @param alpha Alpha.
//' @param lambda A vector representing the intensities of the point processes.
//' @param radius Interaction radius.
//' @param steps Number of steps in the Metropolis algorithm.
//' @param nsim Number of samples to generate.
//' @param types Types of the points.
//' @param model String representing the model to simulate from. At the moment, either "identity", "Strauss" or "Geyer".
//' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
[[nodiscard]] SEXP rmultigibbs(SEXP window, Rcpp::NumericMatrix alpha = 1, Rcpp::NumericVector lambda = 1, double radius = 0, R_xlen_t steps = 30000, R_xlen_t nsim = 1, Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue, Rcpp::CharacterVector model = "identity", bool drop = true) {
  if(Rf_inherits(window, "Rectangle_window")) {
    return rmultigibbs_helper<Window::rectangle>(window, alpha, lambda, radius, steps, nsim, types, model, drop);
  } else if(Rf_inherits(window, "Disk_window")) {
    return rmultigibbs_helper<Window::disk>(window, alpha, lambda, radius, steps, nsim, types, model, drop);
  } else {
    Rcpp::stop("Only rectangle/disk window implemented for now.");
  }
}
