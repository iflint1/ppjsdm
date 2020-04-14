#ifndef INCLUDE_PPJSDM_WINDOW_UTILITIES
#define INCLUDE_PPJSDM_WINDOW_UTILITIES

#include <Rcpp.h>
#include <Rinternals.h>

#include "../point/point_manipulation.hpp"
#include "../utility/im_wrapper.hpp"
#include "window_base.hpp"

#include <cmath> // std::sqrt
#include <string> // std::string

namespace ppjsdm {

class Rectangle_window {
public:
  Rectangle_window(Rcpp::NumericVector marked_range): x_0_(0), delta_x_(1), y_0_(0), delta_y_(1),
  mark_lower_(marked_range[0]),
  delta_mark_(marked_range[1] - marked_range[0]) {}

  explicit Rectangle_window(Rcpp::List window, Rcpp::NumericVector marked_range):
  mark_lower_(marked_range[0]),
  delta_mark_(marked_range[1] - marked_range[0]) {
    const auto x(Rcpp::as<Rcpp::NumericVector>(window["x_range"]));
    x_0_ = x[0];
    delta_x_ = x[1] - x_0_;

    const auto y(Rcpp::as<Rcpp::NumericVector>(window["y_range"]));
    y_0_ = y[0];
    delta_y_ = y[1] - y_0_;
  }

  Marked_point sample(int type = 0) const {
    // TODO: Also sample mark?
    return Marked_point(x_0_ + unif_rand() * delta_x_,  y_0_ + unif_rand() * delta_y_, type, mark_lower_ + delta_mark_ * unif_rand());
  }

  double volume() const {
    return delta_x_ * delta_y_;
  }

  double square_diameter() const {
    return delta_x_ * delta_x_ + delta_y_ * delta_y_;
  }

  double diameter() const {
    return std::sqrt(square_diameter());
  }

private:
  double x_0_;
  double delta_x_;
  double y_0_;
  double delta_y_;
  double mark_lower_;
  double delta_mark_;
};

class Disk_window {
public:
  Disk_window(Rcpp::NumericVector marked_range): x_(0), y_(0), radius_(1),
  mark_lower_(marked_range[0]),
  delta_mark_(marked_range[1] - marked_range[0]) {}

  explicit Disk_window(Rcpp::List window, Rcpp::NumericVector marked_range):
  mark_lower_(marked_range[0]),
  delta_mark_(marked_range[1] - marked_range[0])  {
    const auto centre(Rcpp::as<Rcpp::NumericVector>(window["centre"]));
    x_ = centre[0];
    y_ = centre[1];

    radius_ = static_cast<double>(window["radius"]);
  }

  Marked_point sample(int type = 0) const {
    while(true) {
      const auto x(2 * unif_rand() - 1.0);
      const auto y(2 * unif_rand() - 1.0);
      if(x * x + y * y <= 1) {
        return Marked_point(x_ + radius_ * x,  y_ + radius_ * y, type, mark_lower_ + delta_mark_ * unif_rand());
      }
    }
  }

  double volume() const {
    return M_PI * radius_ * radius_;
  }

  double square_diameter() const {
    return 4 * radius_ * radius_;
  }

  double diameter() const {
    return 2 * radius_;
  }

private:
  double x_;
  double y_;
  double radius_;
  double mark_lower_;
  double delta_mark_;
};

class Im_window {
public:
  explicit Im_window(Rcpp::List im, Rcpp::NumericVector marked_range):
    im_(im),
    volume_(im_.volume()),
    mark_lower_(marked_range[0]),
    delta_mark_(marked_range[1] - marked_range[0]) {}

  Marked_point sample(int type = 0) const {
    const auto x_min(im_.x_min());
    const auto y_min(im_.y_min());
    const auto delta_x(im_.delta_x());
    const auto delta_y(im_.delta_y());
    while(true) {
      const auto x(delta_x * unif_rand() + x_min);
      const auto y(delta_y * unif_rand() + y_min);
      if(im_.is_in(x, y)) {
        return Marked_point(x,  y, type, mark_lower_ + delta_mark_ * unif_rand());
      }
    }
  }

  double volume() const {
    return volume_;
  }

  double square_diameter() const {
    return im_.square_diameter();
  }

  double diameter() const {
    return std::sqrt(square_diameter());
  }

private:
  Im_wrapper im_;
  double volume_;
  double mark_lower_;
  double delta_mark_;
};

inline Window get_window_from_R_object(SEXP window, Rcpp::NumericVector marked_range) {
  if(Rf_isNull(window)) {
    return Window(Rectangle_window(marked_range));
  }
  else {
    const std::string window_class = Rcpp::as<Rcpp::RObject>(window).attr("class");
    if(window_class == "Rectangle_window") {
      return Window(Rectangle_window(window, marked_range));
    } else if(window_class == "Disk_window") {
      return Window(Disk_window(window, marked_range));
    } else if(window_class == "im") {
      return Window(Im_window(window, marked_range));
    } else {
      Rcpp::stop("Unrecognised window type.");
    }
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_WINDOW_UTILITIES
