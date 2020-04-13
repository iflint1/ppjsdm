#ifndef INCLUDE_POTENTIAL
#define INCLUDE_POTENTIAL

#include "../../point/point_manipulation.hpp"
#include "../../point/square_distance.hpp"

namespace ppjsdm {

class Potential {
public:
  static const bool is_nonincreasing;
  static const bool is_nonincreasing_after_lower_endpoint;
  static const bool is_two_valued;
  virtual double apply(double normalized_square_distance, int i, int j) const = 0;
};

template<typename Point, typename Other>
inline auto apply_potential(const Potential& potential, const Point& point, const Other& other) {
  return potential.apply(normalized_square_distance(point, other), get_type(point), get_type(other));
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SHORT_RANGE_POTENTIALS
