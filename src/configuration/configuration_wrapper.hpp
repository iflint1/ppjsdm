#ifndef INCLUDE_PPJSDM_CONFIGURATION_WRAPPER
#define INCLUDE_PPJSDM_CONFIGURATION_WRAPPER

#include <Rcpp.h>
#include <Rinternals.h>

#include "configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"

#include <cstddef> // std::size_t
#include <iterator> // std::iterator

namespace ppjsdm {

class Configuration_wrapper {
private:
  class Marked_point_reference {
  private:
    friend struct traits::point_manipulation<Marked_point_reference>;
    std::tuple<double&, double&, int&> point_;
  public:
    Marked_point_reference(double& x, double& y, int& type):
    point_(x, y, type) {}

    Marked_point_reference& operator=(const Marked_point& point) {
      std::get<0>(point_) = get_x(point);
      std::get<1>(point_) = get_y(point);
      std::get<2>(point_) = get_type(point) + 1;
      return *this;
    }

    operator Marked_point() const {
      return Marked_point(std::get<0>(point_), std::get<1>(point_), std::get<2>(point_) - 1);
    }
  };

  class Const_marked_point_reference {
  private:
    friend struct traits::point_manipulation<Const_marked_point_reference>;
    std::tuple<const double&, const double&, const int&> point_;
  public:
    Const_marked_point_reference(const double& x, const double& y, const int& type):
    point_(x, y, type) {}

    operator Marked_point() const {
      return Marked_point(std::get<0>(point_), std::get<1>(point_), std::get<2>(point_) - 1);
    }
  };

  friend struct traits::point_manipulation<Marked_point_reference>;
  friend struct traits::point_manipulation<Const_marked_point_reference>;
public:
  // TODO: A lot of this iterator and reference stuff is untested; use at your own risk or do some more testing!
  template<typename N, typename I, typename Point>
  class Configuration_iterator : public std::iterator<std::random_access_iterator_tag, Point> {
    friend class Configuration_wrapper;
  protected:
    using Base = std::iterator<std::random_access_iterator_tag, Point>;
    N* px_;
    N* py_;
    I* ptypes_;
    std::size_t i_;

    Configuration_iterator(N* px, N* py, I* ptypes, std::size_t i) :
      px_(px), py_(py), ptypes_(ptypes), i_(i) {}
  public:
    using pointer = typename Base::pointer;
    using reference = std::conditional_t<std::is_const<Point>::value, Const_marked_point_reference, Marked_point_reference>;
    using difference_type = typename Base::difference_type;

    Configuration_iterator(const Configuration_iterator& other) : px_(other.px_), py_(other.py_), ptypes_(other.ptypes_), i_(other.i_) {}

    Configuration_iterator& operator=(const Configuration_iterator& other) {
      px_ = other.px_;
      py_ = other.py_;
      ptypes_ = other.ptypes_;
      i_ = other.i_;
      return *this;
    }

    reference operator*() const {
      return reference((*px_)[i_], (*py_)[i_], (*ptypes_)[i_]);
    }

    Configuration_iterator& operator++() {
      ++i_;
      return *this;
    }

    Configuration_iterator& operator--() {
      --i_;
      return *this;
    }

    Configuration_iterator operator++(int) {
      return Configuration_iterator(px_, py_, ptypes_, i_++);
    }

    Configuration_iterator operator--(int) {
      return Configuration_iterator(px_, py_, ptypes_, i_--);
    }

    Configuration_iterator operator+(difference_type n) const {
      return Configuration_iterator(px_, py_, ptypes_, i_ + n);
    }

    // TODO: Plenty of stuff to factorise here.
    Configuration_iterator& operator+=(difference_type n) {
      i_ += n;
      return *this;
    }

    Configuration_iterator operator-(difference_type n) const {
      return Configuration_iterator(px_, py_, ptypes_, i_ - n);
    }

    Configuration_iterator& operator-=(difference_type n) {
      i_ -= n;
      return *this;
    }

    reference operator[](difference_type n) const {
      return Point((*px_)[i_], (*py_)[i_], (*ptypes_)[i_]);
    }

    bool operator==(const Configuration_iterator& other) const {
      return i_ == other.i_;
    }

    bool operator!=(const Configuration_iterator& other) const {
      return i_ != other.i_;
    }

    bool operator<(const Configuration_iterator& other) const {
      return i_ < other.i_;
    }

    bool operator>(const Configuration_iterator& other) const {
      return i_ > other.i_;
    }

    bool operator<=(const Configuration_iterator& other) const {
      return i_ <= other.i_;
    }

    bool operator>=(const Configuration_iterator& other) const {
      return i_ >= other.i_;
    }

    difference_type operator+(const Configuration_iterator& other) const {
      return i_ + other.i_;
    }

    difference_type operator-(const Configuration_iterator& other) const {
      return i_ - other.i_;
    }
  };

  using iterator = Configuration_iterator<Rcpp::NumericVector, Rcpp::IntegerVector, Marked_point>;
  using const_iterator = Configuration_iterator<const Rcpp::NumericVector, const Rcpp::IntegerVector, const Marked_point>;

  explicit Configuration_wrapper(Rcpp::List configuration):
    x_(Rcpp::as<Rcpp::NumericVector>(configuration["x"])),
    y_(Rcpp::as<Rcpp::NumericVector>(configuration["y"])),
    types_(Rcpp::as<Rcpp::IntegerVector>(configuration["types"])) {}

  explicit Configuration_wrapper(R_xlen_t size):
    x_(Rcpp::no_init(size)),
    y_(Rcpp::no_init(size)),
    types_(Rcpp::no_init(size)) {}
  Configuration_wrapper(): Configuration_wrapper(static_cast<R_xlen_t>(0)) {}

  auto operator[](R_xlen_t index) const {
    return Marked_point(x_[index], y_[index], types_[index] - 1);
  }

  void erase(iterator it) {
    x_.erase(x_[it.i_]);
    y_.erase(y_[it.i_]);
    types_.erase(types_[it.i_]);
  }

  Rcpp::NumericVector x() const {
    return x_;
  }

  Rcpp::NumericVector y() const {
    return y_;
  }

  Rcpp::IntegerVector types() const {
    return types_;
  }

  Marked_point_reference operator[](R_xlen_t index) {
    return Marked_point_reference(x_[index], y_[index], types_[index]);
  }

  auto emplace_back(double x, double y, int type) {
    x_.push_back(x);
    y_.push_back(y);
    types_.push_back(type + 1);
  }

  auto push_back(const Marked_point& point) {
    emplace_back(get_x(point), get_y(point), get_type(point));
  }

  auto size() const {
    return x_.size();
  }

  bool empty() const {
    return size() == 0;
  }

  iterator begin() {
    return iterator(&x_, &y_, &types_, 0);
  }

  const_iterator begin() const {
    return const_iterator(&x_, &y_, &types_, 0);
  }

  const_iterator cbegin() const {
    return const_iterator(&x_, &y_, &types_, 0);
  }

  iterator end() {
    return iterator(&x_, &y_, &types_, x_.size());
  }

  const_iterator end() const {
    return const_iterator(&x_, &y_, &types_, x_.size());
  }

  const_iterator cend() const {
    return const_iterator(&x_, &y_, &types_, x_.size());
  }
private:
  Rcpp::NumericVector x_;
  Rcpp::NumericVector y_;
  Rcpp::IntegerVector types_;
};

namespace traits {

template<>
struct point_manipulation<Configuration_wrapper::Const_marked_point_reference> {
  static inline auto get_x(const Configuration_wrapper::Const_marked_point_reference& point) {
    return std::get<0>(point.point_);
  }
  static inline auto get_y(const Configuration_wrapper::Const_marked_point_reference& point) {
    return std::get<1>(point.point_);
  }
  static inline auto get_type(const Configuration_wrapper::Const_marked_point_reference& point) {
    return std::get<2>(point.point_) - 1;
  }
};

template<>
struct point_manipulation<Configuration_wrapper::Marked_point_reference> {
  static inline auto get_x(const Configuration_wrapper::Marked_point_reference& point) {
    return std::get<0>(point.point_);
  }
  static inline auto get_y(const Configuration_wrapper::Marked_point_reference& point) {
    return std::get<1>(point.point_);
  }
  static inline auto get_type(const Configuration_wrapper::Marked_point_reference& point) {
    return std::get<2>(point.point_) - 1;
  }
};

} // namespace traits
} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_CONFIGURATION_WRAPPER
