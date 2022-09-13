#ifndef INCLUDE_PPJSDM_LIGHTWEIGHT_MATRIX
#define INCLUDE_PPJSDM_LIGHTWEIGHT_MATRIX

#include <Rcpp.h>

#include <type_traits> // std::enable_if, std::is_same
#include <vector> // std::vector

namespace ppjsdm {

template<typename T>
class Lightweight_matrix {
private:
  using MatrixType = std::vector<T>;
public:
  using size_type = typename MatrixType::size_type;
  using value_type = T;

  explicit Lightweight_matrix(size_type rows, size_type columns):
  rows_(rows), columns_(columns), matrix_(rows * columns) {}

  explicit Lightweight_matrix():
    Lightweight_matrix(static_cast<size_type>(0), static_cast<size_type>(0)) {}

  template<typename Matrix, std::enable_if_t<std::is_same<Matrix, Rcpp::NumericMatrix>::value || std::is_same<Matrix, Rcpp::LogicalMatrix>::value || std::is_same<Matrix, Rcpp::IntegerMatrix>::value>* = nullptr>
  explicit Lightweight_matrix(Matrix matrix):
    Lightweight_matrix(matrix.nrow(), matrix.ncol()) {
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

  auto get_row(size_type i) const {
    std::vector<value_type> row(columns_);
    for(size_type j(0); j < columns_; ++j) {
      row[j] = operator()(i, j);
    }
    return row;
  }

  auto get_col(size_type j) const {
    std::vector<value_type> col(rows_);
    for(size_type i(0); i < rows_; ++i) {
      col[i] = operator()(i, j);
    }
    return col;
  }

private:
  size_type rows_;
  size_type columns_;
  MatrixType matrix_;
};

template<typename T>
class Lightweight_square_matrix {
private:
  using MatrixType = std::vector<T>;
public:
  using size_type = typename MatrixType::size_type;
  using value_type = T;

  explicit Lightweight_square_matrix(size_type dim):
    columns_(dim), matrix_(dim * dim) {}

  explicit Lightweight_square_matrix():
    Lightweight_square_matrix(static_cast<size_type>(0)) {}

  template<typename Matrix, std::enable_if_t<std::is_same<Matrix, Rcpp::NumericMatrix>::value || std::is_same<Matrix, Rcpp::LogicalMatrix>::value || std::is_same<Matrix, Rcpp::IntegerMatrix>::value>* = nullptr>
  explicit Lightweight_square_matrix(Matrix matrix):
    Lightweight_square_matrix(matrix.nrow()) {
    if(static_cast<size_type>(matrix.ncol()) != columns_) {
      Rcpp::stop("The matrix is not a square matrix.");
    }
    for(size_type i(0); i < columns_; ++i) {
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
    return columns_;
  }
private:
  size_type columns_;
  MatrixType matrix_;
};

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_LIGHTWEIGHT_MATRIX
