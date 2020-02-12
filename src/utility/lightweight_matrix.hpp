#ifndef INCLUDE_PPJSDM_LIGHTWEIGHT_MATRIX
#define INCLUDE_PPJSDM_LIGHTWEIGHT_MATRIX

#include <Rcpp.h>

#include <vector> // std::vector

namespace ppjsdm {

template<typename T>
class Lightweight_matrix {
private:
  using MatrixType = std::vector<T>;
  using size_t = typename MatrixType::size_type;
public:
  explicit Lightweight_matrix(size_t rows, size_t columns):
  rows_(rows), columns_(columns), matrix_(rows * columns) {}

  explicit Lightweight_matrix(Rcpp::NumericMatrix matrix):
  rows_(matrix.nrow()), columns_(matrix.ncol()), matrix_(rows_ * columns_) {
    for(size_t i(0); i < rows_; ++i) {
      for(size_t j(0); j < columns_; ++j) {
        operator()(i, j) = static_cast<T>(matrix(i, j));
      }
    }
  }

  decltype(auto) operator()(size_t i, size_t j) const {
    return matrix_[i * columns_ + j];
  }

  decltype(auto) operator()(size_t i, size_t j) {
    return matrix_[i * columns_ + j];
  }

  auto ncol() const {
    return columns_;
  }

  auto nrow() const {
    return rows_;
  }
private:
  size_t rows_;
  size_t columns_;
  MatrixType matrix_;
};

// TODO: Factorise with above?
template<typename T>
class Lightweight_square_matrix {
private:
  using MatrixType = std::vector<T>;
  using size_t = typename MatrixType::size_type;
public:
  explicit Lightweight_square_matrix(size_t dim):
  dim_(dim), matrix_(dim * dim) {}

  explicit Lightweight_square_matrix(Rcpp::NumericMatrix matrix):
  dim_(matrix.nrow()), matrix_(dim_ * dim_) {
    if(matrix.ncol() != dim_) {
      Rcpp::stop("The matrix is not a square matrix, as was expected.");
    }
    for(size_t i(0); i < dim_; ++i) {
      for(size_t j(0); j < dim_; ++j) {
        operator()(i, j) = static_cast<T>(matrix(i, j));
      }
    }
  }

  decltype(auto) operator()(size_t i, size_t j) const {
    return matrix_[i * dim_ + j];
  }

  decltype(auto) operator()(size_t i, size_t j) {
    return matrix_[i * dim_ + j];
  }

  auto ncol() const {
    return dim_;
  }

  auto nrow() const {
    return dim_;
  }
private:
  size_t dim_;
  MatrixType matrix_;
};

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_LIGHTWEIGHT_MATRIX
