#ifndef INCLUDE_PPJSDM_GET_LIST_OR_FIRST_ELEMENT
#define INCLUDE_PPJSDM_GET_LIST_OR_FIRST_ELEMENT

#include <Rcpp.h>
#include <Rinternals.h>

namespace ppjsdm {

inline SEXP get_list_or_first_element(Rcpp::List list, bool condition) {
  if(condition) {
    return Rcpp::wrap(list[0]);
  } else {
    return Rcpp::wrap(list);
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_GET_LIST_OR_FIRST_ELEMENT
