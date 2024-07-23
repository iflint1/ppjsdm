// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// compute_papangelou_cpp
Rcpp::NumericVector compute_papangelou_cpp(SEXP configuration, Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::List model, Rcpp::CharacterVector medium_range_model, Rcpp::List alpha, Rcpp::NumericVector beta0, Rcpp::NumericMatrix beta, Rcpp::NumericMatrix gamma, Rcpp::List covariates, Rcpp::List short_range, Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range, R_xlen_t saturation, Rcpp::IntegerVector type, Rcpp::NumericVector mark, int nthreads);
RcppExport SEXP _ppjsdm_compute_papangelou_cpp(SEXP configurationSEXP, SEXP xSEXP, SEXP ySEXP, SEXP modelSEXP, SEXP medium_range_modelSEXP, SEXP alphaSEXP, SEXP beta0SEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP covariatesSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP saturationSEXP, SEXP typeSEXP, SEXP markSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type configuration(configurationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type medium_range_model(medium_range_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type short_range(short_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type medium_range(medium_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type long_range(long_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type type(typeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark(markSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_papangelou_cpp(configuration, x, y, model, medium_range_model, alpha, beta0, beta, gamma, covariates, short_range, medium_range, long_range, saturation, type, mark, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// compute_vcov
Rcpp::List compute_vcov(SEXP configuration, SEXP dummy, SEXP window, Rcpp::List covariates, SEXP model, Rcpp::CharacterVector medium_range_model, Rcpp::List short_range, Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range, R_xlen_t saturation, Rcpp::NumericVector rho, Rcpp::NumericVector theta, Rcpp::NumericMatrix regressors, Rcpp::List data_list, Rcpp::List estimate_alpha, Rcpp::LogicalMatrix estimate_gamma, bool debug, int nthreads, int npoints, bool multiple_windows, std::string dummy_distribution, Rcpp::NumericVector mark_range);
RcppExport SEXP _ppjsdm_compute_vcov(SEXP configurationSEXP, SEXP dummySEXP, SEXP windowSEXP, SEXP covariatesSEXP, SEXP modelSEXP, SEXP medium_range_modelSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP saturationSEXP, SEXP rhoSEXP, SEXP thetaSEXP, SEXP regressorsSEXP, SEXP data_listSEXP, SEXP estimate_alphaSEXP, SEXP estimate_gammaSEXP, SEXP debugSEXP, SEXP nthreadsSEXP, SEXP npointsSEXP, SEXP multiple_windowsSEXP, SEXP dummy_distributionSEXP, SEXP mark_rangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type configuration(configurationSEXP);
    Rcpp::traits::input_parameter< SEXP >::type dummy(dummySEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type medium_range_model(medium_range_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type short_range(short_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type medium_range(medium_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type long_range(long_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type regressors(regressorsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data_list(data_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type estimate_alpha(estimate_alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalMatrix >::type estimate_gamma(estimate_gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< int >::type npoints(npointsSEXP);
    Rcpp::traits::input_parameter< bool >::type multiple_windows(multiple_windowsSEXP);
    Rcpp::traits::input_parameter< std::string >::type dummy_distribution(dummy_distributionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_vcov(configuration, dummy, window, covariates, model, medium_range_model, short_range, medium_range, long_range, saturation, rho, theta, regressors, data_list, estimate_alpha, estimate_gamma, debug, nthreads, npoints, multiple_windows, dummy_distribution, mark_range));
    return rcpp_result_gen;
END_RCPP
}
// compute_S_cpp
Rcpp::NumericMatrix compute_S_cpp(Rcpp::NumericVector rho, Rcpp::NumericVector theta, Rcpp::NumericMatrix regressors, Rcpp::IntegerVector type, int nthreads);
RcppExport SEXP _ppjsdm_compute_S_cpp(SEXP rhoSEXP, SEXP thetaSEXP, SEXP regressorsSEXP, SEXP typeSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type regressors(regressorsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_S_cpp(rho, theta, regressors, type, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// compute_A1_cpp
Rcpp::NumericMatrix compute_A1_cpp(Rcpp::NumericVector rho, Rcpp::NumericVector theta, Rcpp::NumericMatrix regressors, Rcpp::IntegerVector type, int nthreads);
RcppExport SEXP _ppjsdm_compute_A1_cpp(SEXP rhoSEXP, SEXP thetaSEXP, SEXP regressorsSEXP, SEXP typeSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type regressors(regressorsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_A1_cpp(rho, theta, regressors, type, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// compute_A2_plus_A3_cpp
Rcpp::NumericMatrix compute_A2_plus_A3_cpp(SEXP configuration, SEXP window, Rcpp::List covariates, SEXP model, Rcpp::CharacterVector medium_range_model, Rcpp::List short_range, Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range, R_xlen_t saturation, Rcpp::NumericVector rho, Rcpp::NumericVector theta, Rcpp::NumericMatrix regressors, Rcpp::List data_list, Rcpp::List estimate_alpha, Rcpp::LogicalMatrix estimate_gamma, int nthreads, int npoints, bool multiple_windows, Rcpp::NumericVector mark_range, bool debug, int max_executions);
RcppExport SEXP _ppjsdm_compute_A2_plus_A3_cpp(SEXP configurationSEXP, SEXP windowSEXP, SEXP covariatesSEXP, SEXP modelSEXP, SEXP medium_range_modelSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP saturationSEXP, SEXP rhoSEXP, SEXP thetaSEXP, SEXP regressorsSEXP, SEXP data_listSEXP, SEXP estimate_alphaSEXP, SEXP estimate_gammaSEXP, SEXP nthreadsSEXP, SEXP npointsSEXP, SEXP multiple_windowsSEXP, SEXP mark_rangeSEXP, SEXP debugSEXP, SEXP max_executionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type configuration(configurationSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type medium_range_model(medium_range_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type short_range(short_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type medium_range(medium_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type long_range(long_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type regressors(regressorsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data_list(data_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type estimate_alpha(estimate_alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalMatrix >::type estimate_gamma(estimate_gammaSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< int >::type npoints(npointsSEXP);
    Rcpp::traits::input_parameter< bool >::type multiple_windows(multiple_windowsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    Rcpp::traits::input_parameter< int >::type max_executions(max_executionsSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_A2_plus_A3_cpp(configuration, window, covariates, model, medium_range_model, short_range, medium_range, long_range, saturation, rho, theta, regressors, data_list, estimate_alpha, estimate_gamma, nthreads, npoints, multiple_windows, mark_range, debug, max_executions));
    return rcpp_result_gen;
END_RCPP
}
// compute_G2_cpp
Rcpp::NumericMatrix compute_G2_cpp(SEXP configuration, SEXP dummy, SEXP window, Rcpp::List covariates, SEXP model, Rcpp::CharacterVector medium_range_model, Rcpp::List short_range, Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range, Rcpp::IntegerVector type, R_xlen_t saturation, Rcpp::NumericVector rho, Rcpp::NumericVector theta, Rcpp::NumericMatrix regressors, Rcpp::List estimate_alpha, Rcpp::LogicalMatrix estimate_gamma, int nthreads, std::string dummy_distribution, Rcpp::NumericVector mark_range);
RcppExport SEXP _ppjsdm_compute_G2_cpp(SEXP configurationSEXP, SEXP dummySEXP, SEXP windowSEXP, SEXP covariatesSEXP, SEXP modelSEXP, SEXP medium_range_modelSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP typeSEXP, SEXP saturationSEXP, SEXP rhoSEXP, SEXP thetaSEXP, SEXP regressorsSEXP, SEXP estimate_alphaSEXP, SEXP estimate_gammaSEXP, SEXP nthreadsSEXP, SEXP dummy_distributionSEXP, SEXP mark_rangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type configuration(configurationSEXP);
    Rcpp::traits::input_parameter< SEXP >::type dummy(dummySEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type medium_range_model(medium_range_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type short_range(short_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type medium_range(medium_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type long_range(long_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type type(typeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type regressors(regressorsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type estimate_alpha(estimate_alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalMatrix >::type estimate_gamma(estimate_gammaSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< std::string >::type dummy_distribution(dummy_distributionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_G2_cpp(configuration, dummy, window, covariates, model, medium_range_model, short_range, medium_range, long_range, type, saturation, rho, theta, regressors, estimate_alpha, estimate_gamma, nthreads, dummy_distribution, mark_range));
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
// make_types
SEXP make_types(SEXP types, R_xlen_t size, SEXP might_contain_name);
RcppExport SEXP _ppjsdm_make_types(SEXP typesSEXP, SEXP sizeSEXP, SEXP might_contain_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type types(typesSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type might_contain_name(might_contain_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(make_types(types, size, might_contain_name));
    return rcpp_result_gen;
END_RCPP
}
// make_types2
SEXP make_types2(SEXP types, R_xlen_t size, SEXP might_contain_name, SEXP might_contain_name2);
RcppExport SEXP _ppjsdm_make_types2(SEXP typesSEXP, SEXP sizeSEXP, SEXP might_contain_nameSEXP, SEXP might_contain_name2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type types(typesSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type might_contain_name(might_contain_nameSEXP);
    Rcpp::traits::input_parameter< SEXP >::type might_contain_name2(might_contain_name2SEXP);
    rcpp_result_gen = Rcpp::wrap(make_types2(types, size, might_contain_name, might_contain_name2));
    return rcpp_result_gen;
END_RCPP
}
// make_default_types
SEXP make_default_types(R_xlen_t size);
RcppExport SEXP _ppjsdm_make_default_types(SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< R_xlen_t >::type size(sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(make_default_types(size));
    return rcpp_result_gen;
END_RCPP
}
// prepare_gibbsm_data
Rcpp::List prepare_gibbsm_data(Rcpp::List configuration_list, SEXP window, Rcpp::List covariates, Rcpp::List model, Rcpp::CharacterVector medium_range_model, Rcpp::List short_range, SEXP medium_range, SEXP long_range, R_xlen_t saturation, Rcpp::NumericVector mark_range, R_xlen_t max_dummy, R_xlen_t min_dummy, double dummy_factor, Rcpp::List estimate_alpha, Rcpp::LogicalMatrix estimate_gamma, int nthreads, bool debug, std::string dummy_distribution, Rcpp::CharacterVector type_names);
RcppExport SEXP _ppjsdm_prepare_gibbsm_data(SEXP configuration_listSEXP, SEXP windowSEXP, SEXP covariatesSEXP, SEXP modelSEXP, SEXP medium_range_modelSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP saturationSEXP, SEXP mark_rangeSEXP, SEXP max_dummySEXP, SEXP min_dummySEXP, SEXP dummy_factorSEXP, SEXP estimate_alphaSEXP, SEXP estimate_gammaSEXP, SEXP nthreadsSEXP, SEXP debugSEXP, SEXP dummy_distributionSEXP, SEXP type_namesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type configuration_list(configuration_listSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type medium_range_model(medium_range_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type short_range(short_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type medium_range(medium_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type long_range(long_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type max_dummy(max_dummySEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type min_dummy(min_dummySEXP);
    Rcpp::traits::input_parameter< double >::type dummy_factor(dummy_factorSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type estimate_alpha(estimate_alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalMatrix >::type estimate_gamma(estimate_gammaSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    Rcpp::traits::input_parameter< std::string >::type dummy_distribution(dummy_distributionSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type type_names(type_namesSEXP);
    rcpp_result_gen = Rcpp::wrap(prepare_gibbsm_data(configuration_list, window, covariates, model, medium_range_model, short_range, medium_range, long_range, saturation, mark_range, max_dummy, min_dummy, dummy_factor, estimate_alpha, estimate_gamma, nthreads, debug, dummy_distribution, type_names));
    return rcpp_result_gen;
END_RCPP
}
// prepare_gibbsm_data_with_dummy
Rcpp::List prepare_gibbsm_data_with_dummy(Rcpp::List configuration_list, SEXP dummy, SEXP window, Rcpp::List covariates, Rcpp::List model, Rcpp::CharacterVector medium_range_model, Rcpp::List short_range, SEXP medium_range, SEXP long_range, R_xlen_t saturation, Rcpp::NumericVector mark_range, Rcpp::List estimate_alpha, Rcpp::LogicalMatrix estimate_gamma, int nthreads, bool debug, Rcpp::CharacterVector type_names);
RcppExport SEXP _ppjsdm_prepare_gibbsm_data_with_dummy(SEXP configuration_listSEXP, SEXP dummySEXP, SEXP windowSEXP, SEXP covariatesSEXP, SEXP modelSEXP, SEXP medium_range_modelSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP saturationSEXP, SEXP mark_rangeSEXP, SEXP estimate_alphaSEXP, SEXP estimate_gammaSEXP, SEXP nthreadsSEXP, SEXP debugSEXP, SEXP type_namesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type configuration_list(configuration_listSEXP);
    Rcpp::traits::input_parameter< SEXP >::type dummy(dummySEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type medium_range_model(medium_range_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type short_range(short_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type medium_range(medium_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type long_range(long_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type estimate_alpha(estimate_alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalMatrix >::type estimate_gamma(estimate_gammaSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type type_names(type_namesSEXP);
    rcpp_result_gen = Rcpp::wrap(prepare_gibbsm_data_with_dummy(configuration_list, dummy, window, covariates, model, medium_range_model, short_range, medium_range, long_range, saturation, mark_range, estimate_alpha, estimate_gamma, nthreads, debug, type_names));
    return rcpp_result_gen;
END_RCPP
}
// rbinomialpp_cpp
SEXP rbinomialpp_cpp(SEXP window, SEXP n, R_xlen_t nsim, SEXP types, bool drop, Rcpp::NumericVector mark_range);
RcppExport SEXP _ppjsdm_rbinomialpp_cpp(SEXP windowSEXP, SEXP nSEXP, SEXP nsimSEXP, SEXP typesSEXP, SEXP dropSEXP, SEXP mark_rangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type n(nSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< SEXP >::type types(typesSEXP);
    Rcpp::traits::input_parameter< bool >::type drop(dropSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    rcpp_result_gen = Rcpp::wrap(rbinomialpp_cpp(window, n, nsim, types, drop, mark_range));
    return rcpp_result_gen;
END_RCPP
}
// rbirth_cpp
SEXP rbirth_cpp(R_xlen_t nsim, double horizon, Rcpp::CharacterVector types, R_xlen_t nquad, SEXP window, bool drop, R_xlen_t seed, Rcpp::NumericVector mark_range, SEXP starting_configuration, SEXP dummy, int nthreads, SEXP birth_alpha, Rcpp::NumericVector birth_beta0, SEXP birth_covariates, SEXP birth_beta, SEXP birth_gamma, SEXP birth_short_range, SEXP birth_medium_range, SEXP birth_long_range, R_xlen_t birth_saturation, SEXP birth_model, SEXP birth_medium_range_model, SEXP death_alpha, Rcpp::NumericVector death_beta0, SEXP death_covariates, SEXP death_beta, SEXP death_gamma, SEXP death_short_range, SEXP death_medium_range, SEXP death_long_range, R_xlen_t death_saturation, SEXP death_model, SEXP death_medium_range_model);
RcppExport SEXP _ppjsdm_rbirth_cpp(SEXP nsimSEXP, SEXP horizonSEXP, SEXP typesSEXP, SEXP nquadSEXP, SEXP windowSEXP, SEXP dropSEXP, SEXP seedSEXP, SEXP mark_rangeSEXP, SEXP starting_configurationSEXP, SEXP dummySEXP, SEXP nthreadsSEXP, SEXP birth_alphaSEXP, SEXP birth_beta0SEXP, SEXP birth_covariatesSEXP, SEXP birth_betaSEXP, SEXP birth_gammaSEXP, SEXP birth_short_rangeSEXP, SEXP birth_medium_rangeSEXP, SEXP birth_long_rangeSEXP, SEXP birth_saturationSEXP, SEXP birth_modelSEXP, SEXP birth_medium_range_modelSEXP, SEXP death_alphaSEXP, SEXP death_beta0SEXP, SEXP death_covariatesSEXP, SEXP death_betaSEXP, SEXP death_gammaSEXP, SEXP death_short_rangeSEXP, SEXP death_medium_rangeSEXP, SEXP death_long_rangeSEXP, SEXP death_saturationSEXP, SEXP death_modelSEXP, SEXP death_medium_range_modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< R_xlen_t >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< double >::type horizon(horizonSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type types(typesSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type nquad(nquadSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< bool >::type drop(dropSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type starting_configuration(starting_configurationSEXP);
    Rcpp::traits::input_parameter< SEXP >::type dummy(dummySEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type birth_alpha(birth_alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type birth_beta0(birth_beta0SEXP);
    Rcpp::traits::input_parameter< SEXP >::type birth_covariates(birth_covariatesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type birth_beta(birth_betaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type birth_gamma(birth_gammaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type birth_short_range(birth_short_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type birth_medium_range(birth_medium_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type birth_long_range(birth_long_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type birth_saturation(birth_saturationSEXP);
    Rcpp::traits::input_parameter< SEXP >::type birth_model(birth_modelSEXP);
    Rcpp::traits::input_parameter< SEXP >::type birth_medium_range_model(birth_medium_range_modelSEXP);
    Rcpp::traits::input_parameter< SEXP >::type death_alpha(death_alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type death_beta0(death_beta0SEXP);
    Rcpp::traits::input_parameter< SEXP >::type death_covariates(death_covariatesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type death_beta(death_betaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type death_gamma(death_gammaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type death_short_range(death_short_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type death_medium_range(death_medium_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type death_long_range(death_long_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type death_saturation(death_saturationSEXP);
    Rcpp::traits::input_parameter< SEXP >::type death_model(death_modelSEXP);
    Rcpp::traits::input_parameter< SEXP >::type death_medium_range_model(death_medium_range_modelSEXP);
    rcpp_result_gen = Rcpp::wrap(rbirth_cpp(nsim, horizon, types, nquad, window, drop, seed, mark_range, starting_configuration, dummy, nthreads, birth_alpha, birth_beta0, birth_covariates, birth_beta, birth_gamma, birth_short_range, birth_medium_range, birth_long_range, birth_saturation, birth_model, birth_medium_range_model, death_alpha, death_beta0, death_covariates, death_beta, death_gamma, death_short_range, death_medium_range, death_long_range, death_saturation, death_model, death_medium_range_model));
    return rcpp_result_gen;
END_RCPP
}
// rgibbs_cpp
SEXP rgibbs_cpp(SEXP window, SEXP alpha, Rcpp::NumericVector beta0, SEXP covariates, SEXP beta, SEXP gamma, SEXP short_range, SEXP medium_range, SEXP long_range, R_xlen_t saturation, R_xlen_t steps, R_xlen_t nsim, SEXP types, SEXP model, SEXP medium_range_model, bool drop, Rcpp::NumericVector mark_range, Rcpp::IntegerVector only_simulate_these_types, SEXP conditional_configuration, SEXP starting_configuration, R_xlen_t seed, R_xlen_t nthreads, bool debug);
RcppExport SEXP _ppjsdm_rgibbs_cpp(SEXP windowSEXP, SEXP alphaSEXP, SEXP beta0SEXP, SEXP covariatesSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP saturationSEXP, SEXP stepsSEXP, SEXP nsimSEXP, SEXP typesSEXP, SEXP modelSEXP, SEXP medium_range_modelSEXP, SEXP dropSEXP, SEXP mark_rangeSEXP, SEXP only_simulate_these_typesSEXP, SEXP conditional_configurationSEXP, SEXP starting_configurationSEXP, SEXP seedSEXP, SEXP nthreadsSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< SEXP >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type short_range(short_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type medium_range(medium_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type long_range(long_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< SEXP >::type types(typesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type model(modelSEXP);
    Rcpp::traits::input_parameter< SEXP >::type medium_range_model(medium_range_modelSEXP);
    Rcpp::traits::input_parameter< bool >::type drop(dropSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type only_simulate_these_types(only_simulate_these_typesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type conditional_configuration(conditional_configurationSEXP);
    Rcpp::traits::input_parameter< SEXP >::type starting_configuration(starting_configurationSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(rgibbs_cpp(window, alpha, beta0, covariates, beta, gamma, short_range, medium_range, long_range, saturation, steps, nsim, types, model, medium_range_model, drop, mark_range, only_simulate_these_types, conditional_configuration, starting_configuration, seed, nthreads, debug));
    return rcpp_result_gen;
END_RCPP
}
// rppp_cpp
SEXP rppp_cpp(SEXP window, SEXP lambda, R_xlen_t nsim, SEXP types, bool drop, Rcpp::NumericVector mark_range);
RcppExport SEXP _ppjsdm_rppp_cpp(SEXP windowSEXP, SEXP lambdaSEXP, SEXP nsimSEXP, SEXP typesSEXP, SEXP dropSEXP, SEXP mark_rangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< SEXP >::type types(typesSEXP);
    Rcpp::traits::input_parameter< bool >::type drop(dropSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    rcpp_result_gen = Rcpp::wrap(rppp_cpp(window, lambda, nsim, types, drop, mark_range));
    return rcpp_result_gen;
END_RCPP
}
// rstratpp_cpp
SEXP rstratpp_cpp(SEXP window, SEXP delta_x, SEXP delta_y, R_xlen_t nsim, SEXP types, bool drop, Rcpp::NumericVector mark_range);
RcppExport SEXP _ppjsdm_rstratpp_cpp(SEXP windowSEXP, SEXP delta_xSEXP, SEXP delta_ySEXP, SEXP nsimSEXP, SEXP typesSEXP, SEXP dropSEXP, SEXP mark_rangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type delta_y(delta_ySEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< SEXP >::type types(typesSEXP);
    Rcpp::traits::input_parameter< bool >::type drop(dropSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    rcpp_result_gen = Rcpp::wrap(rstratpp_cpp(window, delta_x, delta_y, nsim, types, drop, mark_range));
    return rcpp_result_gen;
END_RCPP
}
// show_short_range_models
Rcpp::CharacterVector show_short_range_models();
RcppExport SEXP _ppjsdm_show_short_range_models() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(show_short_range_models());
    return rcpp_result_gen;
END_RCPP
}
// show_medium_range_models
Rcpp::CharacterVector show_medium_range_models();
RcppExport SEXP _ppjsdm_show_medium_range_models() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(show_medium_range_models());
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP run_testthat_tests(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_ppjsdm_compute_papangelou_cpp", (DL_FUNC) &_ppjsdm_compute_papangelou_cpp, 17},
    {"_ppjsdm_compute_vcov", (DL_FUNC) &_ppjsdm_compute_vcov, 22},
    {"_ppjsdm_compute_S_cpp", (DL_FUNC) &_ppjsdm_compute_S_cpp, 5},
    {"_ppjsdm_compute_A1_cpp", (DL_FUNC) &_ppjsdm_compute_A1_cpp, 5},
    {"_ppjsdm_compute_A2_plus_A3_cpp", (DL_FUNC) &_ppjsdm_compute_A2_plus_A3_cpp, 21},
    {"_ppjsdm_compute_G2_cpp", (DL_FUNC) &_ppjsdm_compute_G2_cpp, 19},
    {"_ppjsdm_has_duplicates", (DL_FUNC) &_ppjsdm_has_duplicates, 1},
    {"_ppjsdm_make_types", (DL_FUNC) &_ppjsdm_make_types, 3},
    {"_ppjsdm_make_types2", (DL_FUNC) &_ppjsdm_make_types2, 4},
    {"_ppjsdm_make_default_types", (DL_FUNC) &_ppjsdm_make_default_types, 1},
    {"_ppjsdm_prepare_gibbsm_data", (DL_FUNC) &_ppjsdm_prepare_gibbsm_data, 19},
    {"_ppjsdm_prepare_gibbsm_data_with_dummy", (DL_FUNC) &_ppjsdm_prepare_gibbsm_data_with_dummy, 16},
    {"_ppjsdm_rbinomialpp_cpp", (DL_FUNC) &_ppjsdm_rbinomialpp_cpp, 6},
    {"_ppjsdm_rbirth_cpp", (DL_FUNC) &_ppjsdm_rbirth_cpp, 33},
    {"_ppjsdm_rgibbs_cpp", (DL_FUNC) &_ppjsdm_rgibbs_cpp, 23},
    {"_ppjsdm_rppp_cpp", (DL_FUNC) &_ppjsdm_rppp_cpp, 6},
    {"_ppjsdm_rstratpp_cpp", (DL_FUNC) &_ppjsdm_rstratpp_cpp, 7},
    {"_ppjsdm_show_short_range_models", (DL_FUNC) &_ppjsdm_show_short_range_models, 0},
    {"_ppjsdm_show_medium_range_models", (DL_FUNC) &_ppjsdm_show_medium_range_models, 0},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_ppjsdm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
