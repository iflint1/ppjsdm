#include <Rcpp.h>
#include <Rmath.h>

#include <cmath> // std::exp
#include <vector> // std::vector

#include "utility/call_on_list_or_vector.h"
#include "utility/compute_phi_dispersion.h"
#include "utility/configuration_manipulation.h"
#include "utility/get_list_or_first_element.h"
#include "utility/make_default_types.h"
#include "utility/make_R_configuration.h"
#include "utility/rbinomialpp_single.h"
#include "utility/remove_random_point.h"
#include "utility/resolve_defaults.h"
#include "utility/window_utilities.h"

template<typename V, typename S>
inline SEXP rmultigibbs_helper(const V& varphi, const S& window, R_xlen_t steps, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, R_xlen_t point_types) {
  const auto volume(window.volume());
  constexpr double prob(0.5);

  Rcpp::List samples(nsim);

  for(R_xlen_t i(0); i < nsim; ++i) {
    // TODO: This can likely be made faster by using one unique vector, since all 3 have the same length.
    // TODO: Make an std_vector_configuration wrapper to make everything more straightforward.
    // TODO: Preallocate with a rough estimate of final size?
    ppjsdm::StdVector_configuration_wrapper points{};

    R_xlen_t total_number(0);

    for(R_xlen_t step(0); step < steps; ++step) {
      const auto u(unif_rand());
      const auto v(unif_rand());
      if(u <= prob) {
        const R_xlen_t type(Rcpp::sample(point_types, 1, false, R_NilValue, false)[0]);
        const auto location_pair(window.sample());
        const Rcpp::NumericVector location{location_pair.first, location_pair.second};

        // TODO: You can avoid taking the exp by reorganising the ratio, and sampling an exponential r.v. instead.
        const auto papangelou(varphi.compute_papangelou(points, location, type, point_types));
        const auto birth_ratio(papangelou * (1 - prob) * volume * point_types / (prob * (1 + total_number)));

        if(v <= birth_ratio) {
          emplace_point(points, location[0], location[1], type + 1);
          ++total_number;
        }
      } else {
        if(total_number != 0) {
          const auto saved_point(remove_random_point(points));
          // TODO: Send the generic point to the function instead.
          const Rcpp::NumericVector saved_location{std::get<0>(saved_point), std::get<1>(saved_point)};
          const R_xlen_t saved_type(std::get<2>(saved_point) - 1);

          const auto papangelou(varphi.compute_papangelou(points, saved_location, saved_type, point_types));
          const auto death_ratio(prob * total_number / (point_types * (1 - prob) * volume * papangelou));
          if(v <= death_ratio) {
            --total_number;
          } else {
            add_point(points, saved_point);
          }
        }
      }
    }

    samples[i] = ppjsdm::make_R_configuration(points, types);
  }

  return ppjsdm::get_list_or_first_element(samples, nsim, drop);
}

//' Sample a multivariate Gibbs point processes
//'
//' @param window Simulation window.
//' @param alpha Repulsion matrix. Default is a square matrix of same size as types, filled with zeroes.
//' @param lambda A vector representing the intensities of the point processes.
//' Default is a vector of same size as types, filled with ones.
//' @param nu A vector representing the dispersion of the number of points.
//' Default is a vector of same size as types, filled with ones.
//' @param radius Interaction radius. Default is zero.
//' @param steps Number of steps in the Metropolis algorithm. Default is 30000.
//' @param nsim Number of samples to generate. Default is 1.
//' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
//' @param model String representing the model to simulate from. At the moment, either "identity", "Strauss", "Geyer" or "neighbour",
//' with default beign "identity".
//' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration. Default is TRUE.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
SEXP rmultigibbs(SEXP window = R_NilValue, SEXP alpha = R_NilValue, SEXP lambda = R_NilValue, SEXP nu = R_NilValue, double radius = 0, R_xlen_t steps = 30000, R_xlen_t nsim = 1, SEXP types = R_NilValue, Rcpp::CharacterVector model = "identity", bool drop = true) {
  const auto point_types(ppjsdm::get_number_types_and_check_conformance(alpha, lambda, nu, types));
  alpha = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(point_types, alpha, 0.);
  lambda = ppjsdm::construct_if_missing<Rcpp::NumericVector>(point_types, lambda, 1.);
  nu = ppjsdm::construct_if_missing<Rcpp::NumericVector>(point_types, nu, 1.);
  types = ppjsdm::make_types(types, point_types, lambda, nu);
  return ppjsdm::call_on_wrapped_window(window, [&alpha, &lambda, &nu, radius, steps, nsim, &types, &model, drop, point_types](const auto& w) {
    return ppjsdm::call_on_list_or_vector(lambda, [&alpha, &lambda, &nu, radius, steps, nsim, &types, &model, drop, point_types, &w](const auto& l) {
      return ppjsdm::call_on_list_or_vector(nu, [&alpha, &lambda, &nu, radius, steps, nsim, &types, &model, drop, point_types, &w, &l](const auto& n) {
        return ppjsdm::call_on_model(model, alpha, l, n, radius, [&w, steps, nsim, &types, drop, point_types](const auto& varphi){
          return rmultigibbs_helper(varphi, w, steps, nsim, types, drop, point_types);
        });
      });
    });
  });
}
