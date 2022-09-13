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
  Window(SEXP window, Rcpp::NumericVector mark_range):
  object_(make_window(window, mark_range)) {}

  explicit Window(SEXP window):
  Window(window, Rcpp::NumericVector::create(1, 1)) {}

  Marked_point sample(std::mt19937 generator, int type) const {
    return object_->sample(generator, type);
  }

  Marked_point sample(int type) const {
    return object_->sample(type);
  }

  double draw_mark() const {
    return object_->draw_mark();
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

  double xmin() const {
    return object_->xmin();
  }

  double xmax() const {
    return object_->xmax();
  }

  double ymin() const {
    return object_->ymin();
  }

  double ymax() const {
    return object_->ymax();
  }

  bool is_in(double x, double y) const {
    return object_->is_in(x, y);
  }

  bool shrink_by_distance(double R) {
    return object_->shrink_by_distance(R);
  }

  bool shrink_by_percent(double percent) {
    return object_->shrink_by_percent(percent);
  }

private:
  struct Concept {
    virtual ~Concept() {}
    // TODO: Do not hard-code generator here and in window_concrete.hpp
    virtual Marked_point sample(std::mt19937 generator, int type) const = 0;
    virtual Marked_point sample(int type) const = 0;
    virtual double draw_mark() const = 0;
    virtual double volume() const = 0;
    virtual double square_diameter() const = 0;
    virtual double diameter() const = 0;
    virtual double xmin() const = 0;
    virtual double xmax() const = 0;
    virtual double ymin() const = 0;
    virtual double ymax() const = 0;
    virtual bool is_in(double x, double y) const = 0;
    virtual bool shrink_by_distance(double R) = 0;
    virtual bool shrink_by_percent(double R) = 0;
  };

  template<typename T>
  class Concrete_window: public Concept {
  public:
    template<typename... Args>
    explicit Concrete_window(Args&&... args): object_(std::forward<Args>(args)...) {}

    Marked_point sample(std::mt19937 generator, int type) const {
      return object_.sample(generator, type);
    }

    Marked_point sample(int type) const {
      return object_.sample(type);
    }

    double draw_mark() const {
      return object_.draw_mark();
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

    double xmin() const {
      return object_.xmin();
    }

    double xmax() const {
      return object_.xmax();
    }

    double ymin() const {
      return object_.ymin();
    }

    double ymax() const {
      return object_.ymax();
    }

    bool is_in(double x, double y) const {
      return object_.is_in(x, y);
    }

    bool shrink_by_distance(double R) {
      return object_.shrink_by_distance(R);
    }

    bool shrink_by_percent(double percent) {
      return object_.shrink_by_percent(percent);
    }

  private:
    T object_;
  };

  static std::shared_ptr<Concept> make_window(SEXP window, Rcpp::NumericVector marked_range) {
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

  std::shared_ptr<Concept> object_;
};

} // namespace ppjsdm

#endif // INCLUDE_WINDOW
