#ifndef INCLUDE_PPJSDM_IS_SYMMETRIC_MATRIX
#define INCLUDE_PPJSDM_IS_SYMMETRIC_MATRIX

#include <Rcpp.h>
#include <Rinternals.h>

namespace ppjsdm {

inline bool is_symmetric_matrix(Rcpp::NumericMatrix matrix) {
  const auto rows(matrix.nrow());
  const auto columns(matrix.ncol());
  for(R_xlen_t i(0); i < rows; ++i) {
    for(R_xlen_t j(i + 1); j < columns; ++j) {
      if(matrix(i, j) != matrix(j, i)) {
        return false;
      }
    }
  }
  return true;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_IS_SYMMETRIC_MATRIX
