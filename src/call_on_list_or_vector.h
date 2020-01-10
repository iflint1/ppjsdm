#ifndef INCLUDE_PPJSDM_CALL_ON_LIST_OR_VECTOR
#define INCLUDE_PPJSDM_CALL_ON_LIST_OR_VECTOR

#include <Rcpp.h>
#include <Rinternals.h>

#include <utility> // std::forward

namespace ppjsdm {

template<typename F>
inline auto call_on_list_or_vector(SEXP generic, F&& f) {
  if(Rf_isNewList(generic)) {
    const Rcpp::List list(generic);
    return std::forward<F>(f)(list);

  } else {
    const Rcpp::NumericVector vector(generic);
    return std::forward<F>(f)(vector);
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CALL_ON_LIST_OR_VECTOR
