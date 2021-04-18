#ifndef INCLUDE_PPJSDM_CALL_ON_LIST_OR_VECTOR
#define INCLUDE_PPJSDM_CALL_ON_LIST_OR_VECTOR

#include <Rcpp.h>
#include <Rinternals.h>

namespace ppjsdm {

template<typename F>
inline auto call_on_list_or_vector(SEXP generic, const F& f) {
  if(Rf_isNewList(generic)) {
    return f(Rcpp::List(generic));
  } else if(Rf_isVector(generic)) {
    return f(Rcpp::NumericVector(generic));
  } else {
    Rcpp::stop("Tried to call call_on_list_or_vector on a SEXP that was neither a list nor a vector.");
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CALL_ON_LIST_OR_VECTOR
