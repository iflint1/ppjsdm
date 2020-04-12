#ifndef INCLUDE_WINDOW_BASE
#define INCLUDE_WINDOW_BASE

namespace ppjsdm {

class Window {
public:
  virtual Marked_point sample(int type) const = 0;
  virtual double volume() const = 0;
  virtual double square_diameter() const = 0;
  virtual double diameter() const = 0;
};

} // namespace ppjsdm

#endif // INCLUDE_WINDOW_BASE
