#ifndef INCLUDE_WINDOW
#define INCLUDE_WINDOW

#include <Rcpp.h>
#include <Rinternals.h>

#include "window_concrete.hpp"

#include "../point/point_manipulation.hpp"

#include <memory> // std::shared_ptr
#include <string> // std::string
#include <utility> // std::forward

namespace ppjsdm {

class Window {
public:
  Window(SEXP window, Rcpp::NumericVector marked_range):
  object_(make_window(window, marked_range)) {}

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
    template<typename... Args>
    explicit Concrete_window(Args&&... args): object_(std::forward<Args>(args)...) {}

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

  static std::shared_ptr<const Concept> make_window(SEXP window, Rcpp::NumericVector marked_range) {
    if(Rf_isNull(window)) {
      return std::make_shared<Concrete_window<detail::Rectangle_window>>(marked_range);
    }
    else {
      const std::string window_class = Rcpp::as<Rcpp::RObject>(window).attr("class");
      if(window_class == "Rectangle_window") {
        return std::make_shared<Concrete_window<detail::Rectangle_window>>(window, marked_range);
      } else if(window_class == "Disk_window") {
        return std::make_shared<Concrete_window<detail::Disk_window>>(window, marked_range);
      } else if(window_class == "im") {
        return std::make_shared<Concrete_window<detail::Im_window>>(window, marked_range);
      } else if(window_class == "Rectangle_window_union") {
        return std::make_shared<Concrete_window<detail::Rectangle_window_union>>(window, marked_range);
      } else {
        Rcpp::stop("Unrecognised window type.");
      }
    }
  }

  std::shared_ptr<const Concept> object_;
};

} // namespace ppjsdm

#endif // INCLUDE_WINDOW
