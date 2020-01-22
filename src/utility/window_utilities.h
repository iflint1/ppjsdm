#ifndef INCLUDE_PPJSDM_WINDOW_UTILITIES
#define INCLUDE_PPJSDM_WINDOW_UTILITIES

#include <Rcpp.h>
#include <Rinternals.h>

#include "../point/point_manipulation.h"
#include "../utility/im_wrapper.h"

#include <string> // std::string

namespace ppjsdm {
namespace detail {

enum class Window {rectangle, disk, im};

template <Window WindowType>
class Window_wrapper;

template <>
class Window_wrapper<Window::rectangle> {
public:
  Window_wrapper(): x_0_(0), delta_x_(1), y_0_(0), delta_y_(1) {}
  explicit Window_wrapper(Rcpp::List window) {
    const auto x(Rcpp::as<Rcpp::NumericVector>(window["x_range"]));
    x_0_ = x[0];
    delta_x_ = x[1] - x_0_;

    const auto y(Rcpp::as<Rcpp::NumericVector>(window["y_range"]));
    y_0_ = y[0];
    delta_y_ = y[1] - y_0_;
  }

  auto sample(int type = 0) const {
    return Marked_point(x_0_ + unif_rand() * delta_x_,  y_0_ + unif_rand() * delta_y_, type);
  }

  double volume() const {
    return delta_x_ * delta_y_;
  }

private:
  double x_0_;
  double delta_x_;
  double y_0_;
  double delta_y_;
};

template <>
class Window_wrapper<Window::disk> {
public:
  Window_wrapper(): x_(0), y_(0), radius_(1) {}
  explicit Window_wrapper(Rcpp::List window) {
    const auto centre(Rcpp::as<Rcpp::NumericVector>(window["centre"]));
    x_ = centre[0];
    y_ = centre[1];

    radius_ = static_cast<double>(window["radius"]);
  }

  auto sample(int type = 0) const {
    while(true) {
      const auto x(2 * unif_rand() - 1.0);
      const auto y(2 * unif_rand() - 1.0);
      if(x * x + y * y <= 1) {
        return Marked_point(x_ + radius_ * x,  y_ + radius_ * y, type);
      }
    }
  }

  double volume() const {
    return M_PI * radius_ * radius_;
  }

private:
  double x_;
  double y_;
  double radius_;
};

template <>
class Window_wrapper<Window::im> {
public:
  explicit Window_wrapper(Rcpp::List im):
    im_(im),
    volume_(im_.volume()) {}

  auto sample(int type = 0) const {
    while(true) {
      const auto x(im_.delta_x() * unif_rand() + im_.x_min());
      const auto y(im_.delta_y() * unif_rand() + im_.y_min());
      if(im_.is_in(x, y)) {
        return Marked_point(x,  y, type);
      }
    }
  }

  double volume() const {
    return volume_;
  }

private:
  Im_wrapper im_;
  double volume_;
};

} // namespace detail

template<typename F>
inline auto call_on_wrapped_window(SEXP window, const F& f) {
  if(Rf_isNull(window)) {
    return f(detail::Window_wrapper<detail::Window::rectangle>{});
  }
  else {
    const std::string window_class = Rcpp::as<Rcpp::RObject>(window).attr("class");
    if(window_class == "Rectangle_window") {
      return f(detail::Window_wrapper<detail::Window::rectangle>(window));
    } else if(window_class == "Disk_window") {
      return f(detail::Window_wrapper<detail::Window::disk>(window));
    } else if(window_class == "im") {
      return f(detail::Window_wrapper<detail::Window::im>(window));
    } else {
      Rcpp::stop("Unrecognised window type.");
    }
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_WINDOW_UTILITIES
