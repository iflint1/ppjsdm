#ifndef INCLUDE_PPJSDM_CALL_ON_LIST_OR_VECTOR
#define INCLUDE_PPJSDM_CALL_ON_LIST_OR_VECTOR

#include <Rcpp.h>
#include <Rinternals.h>

namespace ppjsdm {

template<typename F>
inline auto call_on_list_or_vector(SEXP generic, const F& f) {
  return Rf_isNewList(generic) ?
            f(Rcpp::List(generic)) :
            f(Rcpp::NumericVector(generic));
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CALL_ON_LIST_OR_VECTOR
