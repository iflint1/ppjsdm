#ifndef INCLUDE_PPJSDM_PHI_DISTANCE
#define INCLUDE_PPJSDM_PHI_DISTANCE

#include <cmath> // std::sqrt

enum class Varphi {identity, Strauss, inverse_square};

template<Varphi V>
[[nodiscard]] double compute_phi_distance(double x1, double y1, double x2, double y2, double square_R = 0) {
  const auto delta_x{x2 - x1};
  const auto delta_y{y2 - y1};

  const auto square_distance{delta_x * delta_x + delta_y * delta_y};

  // TODO: Should be if constexpr with C++17
  if(V == Varphi::identity) {
    return std::sqrt(square_distance);
  } else if(V == Varphi::inverse_square) {
    return 1. / square_distance;
  } else {
    if(square_distance <= square_R) {
      return 1.;
    } else {
      return 0.;
    }
  }
}

#endif // INCLUDE_PPJSDM_PHI_DISTANCE
