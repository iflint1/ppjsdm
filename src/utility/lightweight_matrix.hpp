#ifndef INCLUDE_PPJSDM_LIGHTWEIGHT_MATRIX
#define INCLUDE_PPJSDM_LIGHTWEIGHT_MATRIX

#include <Rcpp.h>

#include <vector> // std::vector

namespace ppjsdm {

template<typename T>
class Lightweight_matrix {
private:
  using MatrixType = std::vector<T>;
public:
  using size_type = typename MatrixType::size_type;

  explicit Lightweight_matrix(size_type rows, size_type columns):
  rows_(rows), columns_(columns), matrix_(rows * columns) {}

  explicit Lightweight_matrix(Rcpp::NumericMatrix matrix):
    Lightweight_matrix(matrix.nrow()) {
    for(size_type i(0); i < rows_; ++i) {
      for(size_type j(0); j < columns_; ++j) {
        operator()(i, j) = static_cast<T>(matrix(i, j));
      }
    }
  }

  decltype(auto) operator()(size_type i, size_type j) const {
    return matrix_[i * columns_ + j];
  }

  decltype(auto) operator()(size_type i, size_type j) {
    return matrix_[i * columns_ + j];
  }

  auto ncol() const {
    return columns_;
  }

  auto nrow() const {
    return rows_;
  }
private:
  size_type rows_;
  size_type columns_;
  MatrixType matrix_;
};

// TODO: Factorise with above?
template<typename T>
class Lightweight_square_matrix {
private:
  using MatrixType = std::vector<T>;
public:
  using size_type = typename MatrixType::size_type;

  explicit Lightweight_square_matrix(size_type dim):
  dim_(dim), matrix_(dim * dim) {}

  explicit Lightweight_square_matrix(Rcpp::NumericMatrix matrix):
    Lightweight_square_matrix(matrix.nrow()) {
    if(static_cast<size_type>(matrix.ncol()) != dim_) {
      Rcpp::stop("The matrix is not a square matrix, as was expected.");
    }
    for(size_type i(0); i < dim_; ++i) {
      for(size_type j(0); j < dim_; ++j) {
        operator()(i, j) = static_cast<T>(matrix(i, j));
      }
    }
  }

  decltype(auto) operator()(size_type i, size_type j) const {
    return matrix_[i * dim_ + j];
  }

  decltype(auto) operator()(size_type i, size_type j) {
    return matrix_[i * dim_ + j];
  }

  auto ncol() const {
    return dim_;
  }

  auto nrow() const {
    return dim_;
  }
private:
  size_type dim_;
  MatrixType matrix_;
};

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_LIGHTWEIGHT_MATRIX