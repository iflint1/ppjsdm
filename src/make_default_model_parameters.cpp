#include <Rcpp.h>
#include <Rinternals.h>

#include "utility/construct_if_missing.hpp"
#include "utility/get_number_types.hpp"
#include "utility/is_symmetric_matrix.hpp"
#include "utility/make_default_types.hpp"

// [[Rcpp::export]]
SEXP make_default_model_parameters(Rcpp::List alpha,
                                   SEXP beta0,
                                   SEXP covariates,
                                   SEXP beta,
                                   SEXP gamma,
                                   Rcpp::List short_range,
                                   SEXP medium_range,
                                   SEXP long_range,
                                   SEXP types,
                                   int default_number_types) {
  const auto number_types(ppjsdm::get_number_types_and_check_conformance(default_number_types,
                                                                         alpha,
                                                                         gamma,
                                                                         beta0,
                                                                         short_range,
                                                                         medium_range,
                                                                         long_range,
                                                                         types));
  for(R_xlen_t i(0); i < alpha.size(); ++i) {
    alpha[i] = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(alpha[i], 0., number_types);
    if(!ppjsdm::is_symmetric_matrix(alpha[i])) {
      Rcpp::stop("One of the alphas is not symmetric.");
    }
  }
  gamma = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(gamma, 0., number_types);
  beta0 = ppjsdm::construct_if_missing<Rcpp::NumericVector>(beta0, 0., number_types);
  if(!ppjsdm::is_symmetric_matrix(gamma)) {
    Rcpp::stop("Gamma is not symmetric.");
  }

  const auto beta_nrows(number_types);
  const auto beta_ncols(Rcpp::as<Rcpp::List>(covariates).size());
  beta = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(beta, 1., beta_nrows, beta_ncols);
  if(beta_ncols != 0 && (Rcpp::as<Rcpp::NumericMatrix>(beta).nrow() != beta_nrows
                            || Rcpp::as<Rcpp::NumericMatrix>(beta).ncol() != beta_ncols)) {
    Rcpp::stop("The parameter `beta` does not have the right dimensions.");
  }

  types = ppjsdm::detail::make_types(types, number_types, beta0);
  for(R_xlen_t i(0); i < short_range.size(); ++i) {
    short_range[i] = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(short_range[i], 0.1, number_types);
    if(!ppjsdm::is_symmetric_matrix(short_range[i])) {
      Rcpp::stop("One of the interaction short-range radii matrices is not symmetric.");
    }
  }
  const auto me(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(medium_range, 0., number_types));
  const auto lo(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(long_range, 0., number_types));
  if(!ppjsdm::is_symmetric_matrix(me) || !ppjsdm::is_symmetric_matrix(lo)) {
    Rcpp::stop("One of the interaction radii matrices is not symmetric.");
  }

  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta0") = beta0,
                            Rcpp::Named("covariates") = covariates,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("gamma") = gamma,
                            Rcpp::Named("short_range") = short_range,
                            Rcpp::Named("medium_range") = me,
                            Rcpp::Named("long_range") = lo,
                            Rcpp::Named("types") = types);
}
