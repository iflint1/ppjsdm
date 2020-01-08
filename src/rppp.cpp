#include <Rcpp.h>
#include <Rmath.h>

#include "rbinomialpp_single.h"
#include "make_default_types.h"
#include "window_utilities.h"

template<typename S, typename T>
[[nodiscard]] inline SEXP rppp_helper2(const S& window, T lambda, R_xlen_t nsim, Rcpp::Nullable<Rcpp::CharacterVector> types, bool drop) {
  const auto point_types{lambda.size()};
  const auto volume{window.volume()};

  if(types.isNull()) {
    types = Rcpp::wrap(make_default_types(point_types));
  }

  Rcpp::List samples{nsim};
  auto number_points{Rcpp::IntegerVector(Rcpp::no_init(point_types))};

  for(R_xlen_t i{0}; i < nsim; ++i) {
    R_xlen_t total_number{0};
    for(R_xlen_t j{0}; j < point_types; ++j) {
      const auto points_to_add{R::rpois(volume * static_cast<double>(lambda[j]))};
      number_points[j] = points_to_add;
      total_number += points_to_add;
    }

    samples[i] = rbinomialpp_single(window, number_points, types.as(), point_types, total_number);
  }
  if(nsim == 1 && drop == true) {
    return Rcpp::wrap(samples[0]);
  } else {
    return Rcpp::wrap(samples);
  }
}

template<typename W>
[[nodiscard]] inline SEXP rppp_helper(const W& window, SEXP lambda, R_xlen_t nsim, Rcpp::Nullable<Rcpp::CharacterVector> types, bool drop) {
  if(Rf_isNewList(lambda)) {
    const Rcpp::List lambda_list(lambda);
    return rppp_helper2(window, lambda_list, nsim, types, drop);
  } else {
    const Rcpp::NumericVector lambda_vector(lambda);
    return rppp_helper2(window, lambda_vector, nsim, types, drop);
  }
}

//' Sample a Poisson point processes
//'
//' @param window The window.
//' @param lambda A vector representing the intensities of the multipoint Poisson point processes.
//' @param nsim Number of samples to generate.
//' @param types Types of the points.
//' @param drop If nsim = 1 and drop = TRUE, the result will be a Configuration, rather than a list containing a Configuration.
//' @export
//' @useDynLib ppjsdm
//' @import Rcpp
// [[Rcpp::export]]
[[nodiscard]] SEXP rppp(SEXP window, SEXP lambda = Rcpp::NumericVector::create(1), R_xlen_t nsim = 1, Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue, bool drop = true) {
  return call_on_wrapped_window(window, [&lambda, nsim, &types, drop](const auto& w) { return rppp_helper(w, lambda, nsim, types, drop); });
}
