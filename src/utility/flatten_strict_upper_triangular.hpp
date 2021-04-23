#ifndef INCLUDE_FLATTEN_STRICT_UPPER_TRIANGULAR
#define INCLUDE_FLATTEN_STRICT_UPPER_TRIANGULAR

#include <cmath> // std::floor, std::sqrt
#include <tuple> // std::make_pair

namespace ppjsdm {

// Using linear encoding for the triangular matrices used in the computation of A2/A3 in computation of vcov.
// Reference for the formulas here:
// https://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix
inline auto encode_linear(int i, int j, int n) {
  return (n * (n - 1) - (n - i) * (n - i - 1)) / 2 + j - i - 1;
}

inline auto decode_linear(int k, int n) {
  const auto i(n - 2 - static_cast<int>(std::floor(std::sqrt(-8 * k + 4 * n * (n - 1) - 7) / 2. - 0.5)));
  return std::make_pair(i, k + i + 1 + ((n - i) * (n - i - 1) - n * (n - 1)) / 2);
}

} // namespace ppjsdm

#endif // INCLUDE_FLATTEN_STRICT_UPPER_TRIANGULAR
