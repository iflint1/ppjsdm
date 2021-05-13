#ifndef INCLUDE_POTENTIAL
#define INCLUDE_POTENTIAL

namespace ppjsdm {

template<typename FloatType>
struct Potential {
  virtual ~Potential() {}
  virtual bool is_nonincreasing_after_lower_endpoint() const = 0;
  virtual bool is_two_valued() const = 0;
  virtual FloatType apply(FloatType normalized_square_distance, int i, int j) const = 0;
  virtual FloatType get_square_lower_endpoint(int i, int j) const = 0;
};

} // namespace ppjsdm

#endif // INCLUDE_POTENTIAL
