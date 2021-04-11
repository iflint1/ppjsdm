#ifndef INCLUDE_ALGEBRA
#define INCLUDE_ALGEBRA

#include <type_traits> // std::is_same, std::remove_cv_t, remove_reference_t
#include <vector> // std::vector

namespace ppjsdm {

// Computes (Alpha * beta) [index], for a matrix alpha and a vector beta,
// for all values alpha_{ij} such that f(alpha_{ij}) is true.
// Overload when beta is not a scalar.
template<typename Alpha, typename Beta, typename F, typename std::enable_if_t<!std::is_same<Beta, std::remove_cv_t<std::remove_reference_t<decltype(Alpha{}(0, 0))>>>::value>* = nullptr>
inline auto conditional_matrix_times_vector_at_index(const Alpha& alpha,
                                                     const Beta& beta,
                                                     std::remove_cv_t<decltype(alpha.size())> index,
                                                     F&& f) {
  std::remove_cv_t<std::remove_reference_t<decltype(alpha(0, 0))>> sum(0.);
  const auto ncol(alpha.ncol());
  using size_t = std::remove_cv_t<decltype(ncol)>;
  for(size_t col(0); col < ncol; ++col) {
    const auto b(beta[col]);
    if(b != 0.) {
      const auto a(alpha(index, col));
      if(f(a)) {
        sum += a * b;
      }
    }
  }
  return sum;
}

// Specialization of the above when beta is constant
// Overload for scalar beta.
template<typename Alpha, typename Beta, typename F, typename std::enable_if_t<std::is_same<Beta, std::remove_cv_t<std::remove_reference_t<decltype(Alpha{}(0, 0))>>>::value>* = nullptr>
inline auto conditional_matrix_times_vector_at_index(const Alpha& alpha,
                                                     const Beta& beta,
                                                     std::remove_cv_t<decltype(alpha.size())> index,
                                                     F&& f) {
  std::remove_cv_t<std::remove_reference_t<decltype(alpha(0, 0))>> sum(0.);
  if(beta == 0.) {
    return sum;
  }

  const auto ncol(alpha.ncol());
  using size_t = std::remove_cv_t<decltype(ncol)>;
  for(size_t col(0); col < ncol; ++col) {
    const auto a(alpha(index, col));
    if(f(a)) {
      sum += a;
    }
  }
  // The condition below ensures that 0 * Inf = 0 when beta = Inf and sum = 0.
  if(sum != 0.) {
    return sum * beta;
  } else {
    return 0.;
  }
}

// Computes (Alpha * beta) [index], for a matrix alpha and a vector beta.
template<typename Alpha, typename Beta>
inline auto matrix_times_vector_at_index(const Alpha& alpha,
                                         const Beta& beta,
                                         std::remove_cv_t<decltype(alpha.size())> index) {
  // Removing zero values with a != 0 helps guaramtee that
  // 0 * Inf = 0.
  return conditional_matrix_times_vector_at_index(alpha,
                                                  beta,
                                                  index,
                                                  [](auto a) { return a != 0; });
}

// Computes (Alpha^+ * beta) [index], for a matrix alpha and a vector beta.
template<typename Alpha, typename Beta>
inline auto positive_matrix_times_vector_at_index(const Alpha& alpha,
                                                  const Beta& beta,
                                                  std::remove_cv_t<decltype(alpha.size())> index) {
  return conditional_matrix_times_vector_at_index(alpha,
                                                  beta,
                                                  index,
                                                  [](auto a) { return a > 0; });
}

// Computes (Alpha^- * beta) [index], for a matrix alpha and a vector beta.
template<typename Alpha, typename Beta>
inline auto negative_matrix_times_vector_at_index(const Alpha& alpha,
                                                  const Beta& beta,
                                                  std::remove_cv_t<decltype(alpha.size())> index) {
  return conditional_matrix_times_vector_at_index(alpha,
                                                  beta,
                                                  index,
                                                  [](auto a) { return a < 0; });
}

template<typename Alpha, typename Beta>
inline auto positive_matrix_times_vector(const Alpha& alpha,
                                         const Beta& beta) {
  const auto nrow(alpha.nrow());
  std::vector<std::remove_cv_t<std::remove_reference_t<decltype(alpha(0, 0))>>> result(nrow);
  using size_t = std::remove_cv_t<decltype(nrow)>;
  for(size_t row(0); row < nrow; ++row) {
    result[row] = positive_matrix_times_vector_at_index(alpha, beta, row);
  }
  return result;
}

// Does all elements in the vector Alpha[index, ] satisfy f(alpha_{i,j})?
template<typename Alpha, typename F>
inline bool does_column_satisfy_condition(const Alpha& alpha,
                                          std::remove_cv_t<decltype(alpha.size())> index,
                                          F&& f) {
  const auto ncol(alpha.ncol());
  using size_t = std::remove_cv_t<decltype(ncol)>;
  for(size_t column(0); column < ncol; ++column) {
    if(!f(alpha(index, column))) {
      return false;
    }
  }
  return true;
}

// Is Alpha[index, ] equal to zero?
template<typename Alpha>
inline bool is_column_zero(const Alpha& alpha,
                           std::remove_cv_t<decltype(alpha.size())> index) {
  return does_column_satisfy_condition(alpha, index, [](auto a) { return a == 0.; });
}

// Is Alpha[index, ] nonnegative?
template<typename Alpha>
inline bool is_column_nonnegative(const Alpha& alpha,
                                  std::remove_cv_t<decltype(alpha.size())> index) {
  return does_column_satisfy_condition(alpha, index, [](auto a) { return a >= 0.; });
}

// Is Alpha[index, ] nonpositive?
template<typename Alpha>
inline bool is_column_nonpositive(const Alpha& alpha,
                                  std::remove_cv_t<decltype(alpha.size())> index) {
  return does_column_satisfy_condition(alpha, index, [](auto a) { return a <= 0.; });
}

} // namespace ppjsdm

#endif // INCLUDE_ALGEBRA
