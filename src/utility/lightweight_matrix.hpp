#ifndef INCLUDE_PPJSDM_LIGHTWEIGHT_MATRIX
#define INCLUDE_PPJSDM_LIGHTWEIGHT_MATRIX

#include <vector> // std::vector

namespace ppjsdm {

template<typename T>
class Lightweight_matrix {
private:
  using MatrixType = std::vector<double>;
  using size_t = typename MatrixType::size_type;
public:
  Lightweight_matrix(size_t rows, size_t columns):
  rows_(rows), columns_(columns), matrix_(rows * columns) {}

  decltype(auto) operator()(size_t i, size_t j) const {
    return matrix_[i * static_cast<size_t>(columns_) + j];
  }

  decltype(auto) operator()(size_t i, size_t j) {
    return matrix_[i * static_cast<size_t>(columns_) + j];
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

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_LIGHTWEIGHT_MATRIX
