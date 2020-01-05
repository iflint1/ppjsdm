#ifndef INCLUDE_PPJSDM_WINDOW_UTILITIES
#define INCLUDE_PPJSDM_WINDOW_UTILITIES

#include <Rcpp.h>

#include <tuple> // std::pair

enum class Window {rectangle, disk};

template <Window WindowType>
class Window_wrapper;

template <>
class Window_wrapper<Window::rectangle> {
public:
  explicit Window_wrapper(Rcpp::List window) {
    const auto x{Rcpp::as<Rcpp::NumericVector>(window[0])};
    x_0_ = x[0];
    delta_x_ = x[1] - x_0_;

    const auto y{Rcpp::as<Rcpp::NumericVector>(window[1])};
    y_0_ = y[0];
    delta_y_ = y[1] - y_0_;
  }

  [[nodiscard]] std::pair<double, double> sample() const {
    return std::make_pair(x_0_ + unif_rand() * delta_x_,  y_0_ + unif_rand() * delta_y_);
  }

  [[nodiscard]] double volume() const {
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
  explicit Window_wrapper(Rcpp::List window) {
    const auto centre{Rcpp::as<Rcpp::NumericVector>(window[0])};
    x_ = centre[0];
    y_ = centre[1];

    radius_ = static_cast<double>(window[1]);
  }

  [[nodiscard]] std::pair<double, double> sample() const {
    while(true) {
      const auto x{2 * unif_rand() - 1.0};
      const auto y{2 * unif_rand() - 1.0};
      if(x * x + y * y <= 1) {
        return std::make_pair(x_ + radius_ * x,  y_ + radius_ * y);
      }
    }
  }

  [[nodiscard]] double volume() const {
    return M_PI * radius_ * radius_;
  }

private:
  double x_;
  double y_;
  double radius_;
};

[[nodiscard]] inline Window get_window_type(const Rcpp::List window) {
  if(Rf_inherits(window, "Rectangle_window")) {
    return Window::rectangle;
  } else if(Rf_inherits(window, "Disk_window")) {
    return Window::disk;
  } else {
    Rcpp::stop("Incorrect window type.");
  }
}

#endif // INCLUDE_PPJSDM_WINDOW_UTILITIES
