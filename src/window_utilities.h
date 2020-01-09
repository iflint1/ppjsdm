#ifndef INCLUDE_PPJSDM_WINDOW_UTILITIES
#define INCLUDE_PPJSDM_WINDOW_UTILITIES

#include <Rcpp.h>
#include <Rinternals.h>

#include <tuple> // std::pair
#include <utility> // std::forward

namespace ppjsdm {
namespace detail {

enum class Window {rectangle, disk};

template <Window WindowType>
class Window_wrapper;

template <>
class Window_wrapper<Window::rectangle> {
public:
  Window_wrapper(): x_0_(0), delta_x_(1), y_0_(0), delta_y_(1) {}
  explicit Window_wrapper(Rcpp::List window) {
    const Rcpp::NumericVector x(Rcpp::as<Rcpp::NumericVector>(window["x_range"]));
    x_0_ = x[0];
    delta_x_ = x[1] - x_0_;

    const Rcpp::NumericVector y(Rcpp::as<Rcpp::NumericVector>(window["y_range"]));
    y_0_ = y[0];
    delta_y_ = y[1] - y_0_;
  }

  std::pair<double, double> sample() const {
    return std::make_pair(x_0_ + unif_rand() * delta_x_,  y_0_ + unif_rand() * delta_y_);
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
    const Rcpp::NumericVector centre(Rcpp::as<Rcpp::NumericVector>(window["centre"]));
    x_ = centre[0];
    y_ = centre[1];

    radius_ = static_cast<double>(window["radius"]);
  }

  std::pair<double, double> sample() const {
    while(true) {
      const auto x(2 * unif_rand() - 1.0);
      const auto y(2 * unif_rand() - 1.0);
      if(x * x + y * y <= 1) {
        return std::make_pair(x_ + radius_ * x,  y_ + radius_ * y);
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

} // namespace detail

template<typename F, typename... Args>
inline auto call_on_wrapped_window(SEXP window, F&& f, Args&&... args) {
  if(Rf_isNull(window)) {
    const detail::Window_wrapper<detail::Window::rectangle> wrapped_window{};
    return std::forward<F>(f)(wrapped_window, std::forward<Args>(args)...);
  }
  else if(Rf_inherits(window, "Rectangle_window")) { // TODO: Better to do a switch on attr("class")
    const detail::Window_wrapper<detail::Window::rectangle> wrapped_window(window);
    return std::forward<F>(f)(wrapped_window, std::forward<Args>(args)...);
  } else if(Rf_inherits(window, "Disk_window")) {
    const detail::Window_wrapper<detail::Window::disk> wrapped_window(window);
    return std::forward<F>(f)(wrapped_window, std::forward<Args>(args)...);
  } else {
    Rcpp::stop("Unrecognised window type.");
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_WINDOW_UTILITIES
