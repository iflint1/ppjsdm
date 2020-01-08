#ifndef INCLUDE_PPJSDM_CALL_ON_LIST_OR_VECTOR
#define INCLUDE_PPJSDM_CALL_ON_LIST_OR_VECTOR

#include <Rcpp.h>

#include <utility> // std::forward

namespace ppjsdm {

template<typename F, typename... Args>
inline auto call_on_list_or_vector(SEXP generic, F&& f, Args&&... args) {
  if(Rf_isNewList(generic)) {
    const Rcpp::List list(generic);
    return std::forward<F>(f)(list, std::forward<Args>(args)...);

  } else {
    const Rcpp::NumericVector vector(generic);
    return std::forward<F>(f)(vector, std::forward<Args>(args)...);
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CALL_ON_LIST_OR_VECTOR
