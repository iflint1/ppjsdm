#ifndef INCLUDE_WINDOW_BASE
#define INCLUDE_WINDOW_BASE

#include <memory> // std::shared_ptr
#include <type_traits> // std::remove_cv, std::remove_reference
#include <utility> // std::forward

namespace ppjsdm {

class Window {
public:
  template<typename T>
  Window(T&& object):
    object_(std::make_shared<Concrete_window<std::remove_cv_t<std::remove_reference_t<T>>>>(std::forward<T>(object))) {}

  Marked_point sample(int type) const {
    return object_->sample(type);
  }

  double volume() const {
    return object_->volume();
  }

  double square_diameter() const {
    return object_->square_diameter();
  }

  double diameter() const {
    return object_->diameter();
  }

private:
  struct Concept {
    virtual ~Concept() {}
    virtual Marked_point sample(int type) const = 0;
    virtual double volume() const = 0;
    virtual double square_diameter() const = 0;
    virtual double diameter() const = 0;
  };

  template<typename T>
  class Concrete_window: public Concept {
  public:
    template<typename S>
    explicit Concrete_window(S&& object): object_(std::forward<S>(object)) {}

    Marked_point sample(int type) const {
      return object_.sample(type);
    }

    double volume() const {
      return object_.volume();
    }

    double square_diameter() const {
      return object_.square_diameter();
    }

    double diameter() const {
      return object_.diameter();
    }

  private:
    T object_;
  };

  std::shared_ptr<const Concept> object_;
};
} // namespace ppjsdm

#endif // INCLUDE_WINDOW_BASE
