#ifndef INCLUDE_PPJSDM_WINDOW_CONCRETE
#define INCLUDE_PPJSDM_WINDOW_CONCRETE

#include <Rcpp.h>
#include <Rinternals.h>

#include "../point/point_manipulation.hpp"
#include "../utility/im_wrapper.hpp"

#include <cmath> // std::sqrt
#include <vector> // std::vector

namespace ppjsdm {
namespace detail {

class Rectangle_window {
public:
  explicit Rectangle_window(Rcpp::NumericVector marked_range): x_0_(0), delta_x_(1), y_0_(0), delta_y_(1),
  mark_lower_(marked_range[0]),
  delta_mark_(marked_range[1] - marked_range[0]) {}

  explicit Rectangle_window(Rcpp::List window, Rcpp::NumericVector marked_range):
  mark_lower_(marked_range[0]),
  delta_mark_(marked_range[1] - marked_range[0]) {
    const auto x(Rcpp::as<Rcpp::NumericVector>(window["x_range"]));
    x_0_ = x[0];
    delta_x_ = x[1] - x[0];

    const auto y(Rcpp::as<Rcpp::NumericVector>(window["y_range"]));
    y_0_ = y[0];
    delta_y_ = y[1] - y[0];
  }

  Marked_point sample(int type = 0) const {
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

  // TODO: Test is_in here and in other window classes
  template<typename Point>
  bool is_in(const Point& point) const {
    return get_x(point) >= x_0_ && get_x(point) <= x_0_ + delta_x_ && get_y(point) >= y_0_ && get_x(point) <= y_0_ + delta_y_;
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
  explicit Disk_window(Rcpp::NumericVector marked_range): x_(0), y_(0), radius_(1),
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

  template<typename Point>
  bool is_in(const Point& point) const {
    return (get_x(point) - x_) * (get_x(point) - x_) + (get_y(point) - y_) * (get_y(point) - y_) <= radius_ * radius_;
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

  template<typename Point>
  bool is_in(const Point& point) const {
    return im_.is_in(get_x(point), get_y(point));
  }

private:
  Im_wrapper im_;
  double volume_;
  double mark_lower_;
  double delta_mark_;
};

class Rectangle_window_union {
public:
  explicit Rectangle_window_union(Rcpp::NumericVector marked_range): x_0_{0}, delta_x_{1}, y_0_{0}, delta_y_{1},
  mark_lower_(marked_range[0]),
  delta_mark_(marked_range[1] - marked_range[0]) {}

  explicit Rectangle_window_union(Rcpp::List window, Rcpp::NumericVector marked_range):
    mark_lower_(marked_range[0]),
    delta_mark_(marked_range[1] - marked_range[0]) {
    const auto x_ranges(Rcpp::as<Rcpp::List>(window["x_ranges"]));
    const auto n(x_ranges.size());
    x_0_ = std::vector<double>(n);
    delta_x_ = std::vector<double>(n);
    y_0_ = std::vector<double>(n);
    delta_y_ = std::vector<double>(n);
    const auto y_ranges(Rcpp::as<Rcpp::List>(window["y_ranges"]));
    using size_t = decltype(x_ranges.size());
    for(size_t i(0); i < n; ++i) {
      const auto x(Rcpp::as<Rcpp::NumericVector>(x_ranges[i]));
      x_0_[i] = x[0];
      delta_x_[i] = x[1] - x[0];

      const auto y(Rcpp::as<Rcpp::NumericVector>(y_ranges[i]));
      y_0_[i] = y[0];
      delta_y_[i] = y[1] - y[0];
    }
  }

  double volume() const {
    double volume(0);
    using size_t = decltype(delta_x_.size());
    for(size_t i(0); i < delta_x_.size(); ++i) {
      volume += delta_x_[i] * delta_y_[i];
    }
    return volume;
  }

  Marked_point sample(int type = 0) const {
    double total_weight(volume());
    using size_t = decltype(delta_x_.size());
    size_t index(0);
    const auto u(unif_rand());
    auto sum(delta_x_[index] * delta_y_[index]);
    while(sum < u * total_weight) {
      ++index;
      sum += delta_x_[index] * delta_y_[index];
    }
    return Marked_point(x_0_[index] + unif_rand() * delta_x_[index],
                        y_0_[index] + unif_rand() * delta_y_[index],
                                                            type,
                                                            mark_lower_ + delta_mark_ * unif_rand());
  }

  double square_diameter() const {
    // TODO: I suspect I don't need this function here and in other windows
    return 1.;
  }

  double diameter() const {
    // TODO: I suspect I don't need this function here and in other windows
    return 1.;
  }

  template<typename Point>
  bool is_in(const Point& point) const {
    const auto n(x_0_.size());
    using size_t = decltype(x_0_.size());
    for(size_t i(0); i < n; ++i) {
      if(get_x(point) >= x_0_[i] && get_x(point) <= x_0_[i] + delta_x_[i] && get_y(point) >= y_0_[i] && get_x(point) <= y_0_[i] + delta_y_[i]) {
        return true;
      }
    }
    return false;
  }

private:
  // TODO: Group all 4 vectors into a tuple so that you don't have to assume they have the same size
  std::vector<double> x_0_;
  std::vector<double> delta_x_;
  std::vector<double> y_0_;
  std::vector<double> delta_y_;
  double mark_lower_;
  double delta_mark_;
};

} // namespace detail
} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_WINDOW_CONCRETE
