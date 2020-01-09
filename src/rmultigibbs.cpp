#include <Rcpp.h>
#include <Rmath.h>

#include <cmath> // std::exp
#include <vector> // std::vector

#include "compute_phi_dispersion.h"
#include "configuration_utilities.h"
#include "get_list_or_first_element.h"
#include "rbinomialpp_single.h"
#include "make_default_types.h"
#include "window_utilities.h"

namespace ppjsdm {

template<typename V, typename S>
inline SEXP rmultigibbs_helper2(const S& window, R_xlen_t steps, R_xlen_t nsim, Rcpp::CharacterVector types, bool drop, const V& varphi, R_xlen_t point_types) {
  const auto volume(window.volume());
  const double prob(0.5);

  Rcpp::List samples(nsim);

  for(R_xlen_t i(0); i < nsim; ++i) {
    // TODO: This can likely be made faster by using one unique vector, since all 3 have the same length.
    // TODO: Make an std_vector_configuration wrapper to make everything more straightforward.
    // TODO: Preallocate with a rough estimate of final size?
    std::vector<double> x;
    std::vector<double> y;
    std::vector<int> types_vector;

    R_xlen_t total_number(0);

    for(R_xlen_t step(0); step < steps; ++step) {
      const auto u(unif_rand());
      const auto v(unif_rand());
      if(u <= prob) {
        const R_xlen_t type(Rcpp::sample(point_types, 1, false, R_NilValue, false)[0]);
        const auto location_pair(window.sample());
        const Rcpp::NumericVector location{location_pair.first, location_pair.second};

        // TODO: You can avoid taking the exp by reorganising the ratio, and sampling an exponential r.v. instead.
        const auto papangelou(varphi.compute_papangelou(x, y, types_vector, location, type, point_types));
        const auto birth_ratio(papangelou * (1 - prob) * volume * point_types / (prob * (1 + total_number)));

        if(v <= birth_ratio) {
          // Add point to configuration
          x.push_back(location[0]);
          y.push_back(location[1]);
          types_vector.push_back(type + 1);
          ++total_number;
        }
      } else {
        if(total_number != 0) {
          const R_xlen_t index(Rcpp::sample(x.size(), 1, false, R_NilValue, false)[0]);
          const Rcpp::NumericVector saved_location{x[index], y[index]};
          const R_xlen_t saved_type(types_vector[index] - 1);
          x.erase(x.begin() + index);
          y.erase(y.begin() + index);
          types_vector.erase(types_vector.begin() + index);

          const auto papangelou(varphi.compute_papangelou(x, y, types_vector, saved_location, saved_type, point_types));
          const auto death_ratio(prob * total_number / (point_types * (1 - prob) * volume * papangelou));
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

    samples[i] = make_configuration(Rcpp::wrap(x), Rcpp::wrap(y), Rcpp::wrap(types_vector), types);
  }

  return get_list_or_first_element(samples, nsim, drop);
}

template<typename W>
inline SEXP rmultigibbs_helper(const W& window, Rcpp::NumericMatrix alpha, Rcpp::NumericVector lambda, Rcpp::NumericVector nu, double radius, R_xlen_t steps, R_xlen_t nsim, Rcpp::Nullable<Rcpp::CharacterVector> types, Rcpp::CharacterVector model, bool drop) {
  // TODO: Use dispatcher as in windows_utilities.h?
  const auto point_types(lambda.size());
  const auto types_vector(ppjsdm::make_default_types(types, lambda, point_types));
  if(model[0] == "identity") {
    const Exponential_family_model<Varphi_model_papangelou<varphi::Identity>,
                                   decltype(lambda), decltype(nu), decltype(alpha)> varphi(lambda, nu, alpha);
    return rmultigibbs_helper2(window, steps, nsim, types_vector, drop, varphi, point_types);
  } else if(model[0] == "Strauss") {
    const Exponential_family_model<Varphi_model_papangelou<varphi::Strauss>,
                                   decltype(lambda), decltype(nu), decltype(alpha)> varphi(lambda, nu, alpha, radius);
    return rmultigibbs_helper2(window, steps, nsim, types_vector, drop, varphi, point_types);
  } else if(model[0] == "Geyer") {
    const Exponential_family_model<Geyer_papangelou,
                                   decltype(lambda), decltype(nu), decltype(alpha)> varphi(lambda, nu, alpha, radius, 2.0);
    return rmultigibbs_helper2(window, steps, nsim, types_vector, drop, varphi, point_types);
  } else if(model[0] == "neighbour") {
    const Exponential_family_model<Nearest_neighbour_papangelou<varphi::Identity>,
                                   decltype(lambda), decltype(nu), decltype(alpha)> varphi(lambda, nu, alpha);
    return rmultigibbs_helper2(window, steps, nsim, types_vector, drop, varphi, point_types);
  } else {
    Rcpp::stop("Incorrect model entered.\n");
  }
}

} // namespace ppjsdm

//' Sample a multivariate Gibbs point processes
//'
//' @param window The window.
//' @param alpha Alpha.
//' @param lambda A vector representing the intensities of the point processes.
//' @param nu A vector representing the dispersion of the number of points.
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
SEXP rmultigibbs(SEXP window, Rcpp::NumericMatrix alpha, Rcpp::NumericVector lambda, Rcpp::NumericVector nu, double radius = 0, R_xlen_t steps = 30000, R_xlen_t nsim = 1, Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue, Rcpp::CharacterVector model = "identity", bool drop = true) {
  // TODO: Add check for size of alpha.
  // TODO: I'm getting a crash on multiple point types -- debug.
  return ppjsdm::call_on_wrapped_window(window, [&alpha, &lambda, &nu, radius, steps, nsim, &types, &model, drop](const auto& w) {
    return ppjsdm::rmultigibbs_helper(w, alpha, lambda, nu, radius, steps, nsim, types, model, drop);
  });
}
