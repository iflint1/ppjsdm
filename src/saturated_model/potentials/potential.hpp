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

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SHORT_RANGE_POTENTIALS
