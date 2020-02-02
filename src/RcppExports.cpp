// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// compute_papangelou
double compute_papangelou(SEXP configuration, Rcpp::NumericVector coordinates, R_xlen_t type, Rcpp::CharacterVector model, Rcpp::NumericMatrix alpha, Rcpp::NumericVector lambda, Rcpp::NumericMatrix coefs, Rcpp::List covariates, Rcpp::NumericMatrix radius, R_xlen_t saturation);
RcppExport SEXP _ppjsdm_compute_papangelou(SEXP configurationSEXP, SEXP coordinatesSEXP, SEXP typeSEXP, SEXP modelSEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP coefsSEXP, SEXP covariatesSEXP, SEXP radiusSEXP, SEXP saturationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type configuration(configurationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type coordinates(coordinatesSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type type(typeSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type coefs(coefsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_papangelou(configuration, coordinates, type, model, alpha, lambda, coefs, covariates, radius, saturation));
    return rcpp_result_gen;
END_RCPP
}
// has_duplicates
bool has_duplicates(Rcpp::List configuration);
RcppExport SEXP _ppjsdm_has_duplicates(SEXP configurationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type configuration(configurationSEXP);
    rcpp_result_gen = Rcpp::wrap(has_duplicates(configuration));
    return rcpp_result_gen;
END_RCPP
}
// prepare_gibbsm_data
Rcpp::List prepare_gibbsm_data(SEXP configuration, SEXP window, Rcpp::List covariates, Rcpp::CharacterVector model, SEXP radius, R_xlen_t saturation);
RcppExport SEXP _ppjsdm_prepare_gibbsm_data(SEXP configurationSEXP, SEXP windowSEXP, SEXP covariatesSEXP, SEXP modelSEXP, SEXP radiusSEXP, SEXP saturationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type configuration(configurationSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< SEXP >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    rcpp_result_gen = Rcpp::wrap(prepare_gibbsm_data(configuration, window, covariates, model, radius, saturation));
    return rcpp_result_gen;
END_RCPP
}
// rbinomialpp
SEXP rbinomialpp(SEXP window, SEXP n, R_xlen_t nsim, SEXP types, bool drop);
RcppExport SEXP _ppjsdm_rbinomialpp(SEXP windowSEXP, SEXP nSEXP, SEXP nsimSEXP, SEXP typesSEXP, SEXP dropSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type n(nSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< SEXP >::type types(typesSEXP);
    Rcpp::traits::input_parameter< bool >::type drop(dropSEXP);
    rcpp_result_gen = Rcpp::wrap(rbinomialpp(window, n, nsim, types, drop));
    return rcpp_result_gen;
END_RCPP
}
// rgibbs_cpp
SEXP rgibbs_cpp(SEXP window, SEXP alpha, SEXP lambda, SEXP covariates, SEXP coefs, SEXP radius, R_xlen_t saturation, R_xlen_t steps, R_xlen_t nsim, SEXP types, Rcpp::CharacterVector model, bool drop);
RcppExport SEXP _ppjsdm_rgibbs_cpp(SEXP windowSEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP covariatesSEXP, SEXP coefsSEXP, SEXP radiusSEXP, SEXP saturationSEXP, SEXP stepsSEXP, SEXP nsimSEXP, SEXP typesSEXP, SEXP modelSEXP, SEXP dropSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type coefs(coefsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< SEXP >::type types(typesSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< bool >::type drop(dropSEXP);
    rcpp_result_gen = Rcpp::wrap(rgibbs_cpp(window, alpha, lambda, covariates, coefs, radius, saturation, steps, nsim, types, model, drop));
    return rcpp_result_gen;
END_RCPP
}
// rppp
SEXP rppp(SEXP window, SEXP lambda, R_xlen_t nsim, SEXP types, bool drop);
RcppExport SEXP _ppjsdm_rppp(SEXP windowSEXP, SEXP lambdaSEXP, SEXP nsimSEXP, SEXP typesSEXP, SEXP dropSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< SEXP >::type types(typesSEXP);
    Rcpp::traits::input_parameter< bool >::type drop(dropSEXP);
    rcpp_result_gen = Rcpp::wrap(rppp(window, lambda, nsim, types, drop));
    return rcpp_result_gen;
END_RCPP
}
// show_models
void show_models();
RcppExport SEXP _ppjsdm_show_models() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    show_models();
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ppjsdm_compute_papangelou", (DL_FUNC) &_ppjsdm_compute_papangelou, 10},
    {"_ppjsdm_has_duplicates", (DL_FUNC) &_ppjsdm_has_duplicates, 1},
    {"_ppjsdm_prepare_gibbsm_data", (DL_FUNC) &_ppjsdm_prepare_gibbsm_data, 6},
    {"_ppjsdm_rbinomialpp", (DL_FUNC) &_ppjsdm_rbinomialpp, 5},
    {"_ppjsdm_rgibbs_cpp", (DL_FUNC) &_ppjsdm_rgibbs_cpp, 12},
    {"_ppjsdm_rppp", (DL_FUNC) &_ppjsdm_rppp, 5},
    {"_ppjsdm_show_models", (DL_FUNC) &_ppjsdm_show_models, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_ppjsdm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
