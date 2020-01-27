#ifndef INCLUDE_IM_WRAPPER
#define INCLUDE_IM_WRAPPER

#include <Rcpp.h>
#include <Rinternals.h>

#include "../utility/size_t.hpp"

#include <algorithm> // std::max, std::min
#include <cmath> // std::round
#include <limits> // std::numeric_limits
#include <string> // std::string
#include <tuple> // std::pair
#include <vector> // std::vector

namespace ppjsdm {

class Im_wrapper {
private:
  void set_matrix(R_xlen_t i, R_xlen_t j, double value) {
    mat_[i * number_col_ + j] = value;
  }
  auto get_matrix(R_xlen_t i, R_xlen_t j) const {
    return mat_[i * number_col_ + j];
  }
public:
  explicit Im_wrapper(Rcpp::List im):
    number_row_(Rcpp::as<Rcpp::IntegerVector>(im["dim"])[0]),
    number_col_(Rcpp::as<Rcpp::IntegerVector>(im["dim"])[1]),
    xcol_(Rf_asReal(im["xcol"])),
    yrow_(Rf_asReal(im["yrow"])),
    xstep_(Rf_asReal(im["xstep"])),
    ystep_(Rf_asReal(im["ystep"])),
    mat_(number_row_ * number_col_) {
    if(!Rf_inherits(im, "im")) {
      Rcpp::stop("Tried to construct Im_wrapper with an Rcpp::List which is not an im.");
    }
    const auto mat(Rcpp::as<Rcpp::NumericMatrix>(im["v"]));
    for(R_xlen_t i(0); i < number_row_; ++i) {
      for(R_xlen_t j(0); j < number_col_; ++j) {
        set_matrix(i, j, mat(i, j));
      }
    }
  }

  double operator()(double x, double y) const {
    const R_xlen_t raw_col(std::round((x - xcol_) / xstep_));
    const R_xlen_t raw_row(std::round((y - yrow_) / ystep_));
    const R_xlen_t col(std::max<int>(0, std::min<int>(raw_col, number_col_ - 1)));
    const R_xlen_t row(std::max<int>(0, std::min<int>(raw_row, number_row_ - 1)));
    return get_matrix(row, col);
  }

  Rcpp::NumericVector operator()(Rcpp::NumericVector x, Rcpp::NumericVector y) const {
    const R_xlen_t length_x(x.size());
    Rcpp::NumericVector result(Rcpp::no_init(length_x));

    for(R_xlen_t i(0); i < length_x; ++i) {
      result[i] = operator()(x[i], y[i]);
    }
    return result;
  }



  bool is_in(double x, double y) const {
    const auto mat_value(operator()(x, y));
    return !R_IsNA(mat_value);
  }

  double volume() const {
    R_xlen_t non_na_values(0);
    for(R_xlen_t i(0); i < number_row_; ++i) {
      for(R_xlen_t j(0); j < number_col_; ++j) {
        if(!R_IsNA(get_matrix(i, j))) {
          ++non_na_values;
        }
      }
    }
    return static_cast<double>(non_na_values) * xstep_ * ystep_;
  }

  std::pair<double, double> bounds() const {
    double lower(std::numeric_limits<double>::infinity());
    double upper(-std::numeric_limits<double>::infinity());
    for(R_xlen_t i(0); i < number_row_; ++i) {
      for(R_xlen_t j(0); j < number_col_; ++j) {
        const auto value(get_matrix(i, j));
        if(!R_IsNA(value)) {
          if(value < lower) {
            lower = value;
          }
          if(value > lower) {
            upper = value;
          }
        }
      }
    }
    return std::make_pair(lower, upper);
  }

  double diameter() const {
    double square_diameter(0);
    for(R_xlen_t i(0); i < number_row_; ++i) {
      for(R_xlen_t j(0); j < number_col_; ++j) {
        if(!R_IsNA(get_matrix(i, j))) {
          for(R_xlen_t k(0); k < number_row_; ++k) {
            for(R_xlen_t l(0); l < number_col_; ++l) {
              if(!R_IsNA(get_matrix(k, l))) {
                const auto delta_x((i - k) * xstep_);
                const auto delta_y((j - l) * ystep_);
                const auto square_distance(delta_x * delta_x + delta_y * delta_y);
                if(square_distance > square_diameter) {
                  square_diameter = square_distance;
                }
              }
            }
          }
        }
      }
    }
    return std::sqrt(square_diameter);
  }

  double x_min() const {
    return xcol_;
  }

  double delta_x() const {
    return static_cast<double>(number_col_) * xstep_;
  }

  double y_min() const {
    return yrow_;
  }

  double delta_y() const {
    return static_cast<double>(number_row_) * ystep_;
  }
private:
  R_xlen_t number_row_;
  R_xlen_t number_col_;
  double xcol_;
  double yrow_;
  double xstep_;
  double ystep_;
  std::vector<double> mat_;
};

class Im_list_wrapper {
public:
  // Note: Not preallocating im_list_ with im_list_(im_list.size) since each element
  // of the vector will have variable size.
  // Instead, emplace_back when size is known.
  explicit Im_list_wrapper(Rcpp::List im_list) {
    const auto n(im_list.size());
    using size_t = size_t<Rcpp::List>;
    if(n > 0) {
      const auto names = Rcpp::as<Rcpp::CharacterVector>(im_list.names());
      for(size_t i(0); i < n; ++i) {
        im_list_.emplace_back(Rcpp::as<Rcpp::List>(im_list[i]));
        im_names_.emplace_back(names[i]);
      }
    }
  }

  decltype(auto) operator[](R_xlen_t index) {
    return im_list_[index];
  }

  decltype(auto) operator[](R_xlen_t index) const {
    return im_list_[index];
  }

  R_xlen_t size() const {
    return im_list_.size();
  }

  const auto& names() const {
    return im_names_;
  }
private:
  std::vector<Im_wrapper> im_list_;
  std::vector<std::string> im_names_;
};

} // namespace ppjsdm

#endif // INCLUDE_IM_WRAPPER