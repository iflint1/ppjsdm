#ifndef INCLUDE_SHORT_RANGE_POTENTIALS
#define INCLUDE_SHORT_RANGE_POTENTIALS

#include <Rcpp.h>

#include "potential.hpp"
#include "../../point/point_manipulation.hpp"
#include "../../point/square_distance.hpp"
#include "../../utility/lightweight_matrix.hpp"

#include <algorithm> // std::max
#include <cmath> // std::sqrt, std::exp, std::log
#include <memory> // std::shared_ptr

namespace ppjsdm {
namespace potentials {

// Note: Public inheritance in order to inherit member variables if they exist in Varphi.
template<typename Varphi>
class Short_range_potential: public Varphi, public Potential {
private:
  Lightweight_square_matrix<double> matrix_;
  using size_t = typename decltype(matrix_)::size_type;
public:
  explicit Short_range_potential(Rcpp::NumericMatrix radius): matrix_(radius.nrow()) {
    const size_t dim(radius.nrow());
    if(static_cast<size_t>(radius.ncol()) != dim) {
      Rcpp::stop("The matrix is not a square matrix, as was expected.");
    }
    for(size_t i(0); i < static_cast<size_t>(dim); ++i) {
      for(size_t j(0); j < static_cast<size_t>(dim); ++j) {
        matrix_(i, j) = Varphi::set(radius(i, j));
      }
    }
  }

  double apply(double normalized_square_distance, int i, int j) const override {
    return Varphi::apply(normalized_square_distance, matrix_(i, j));
  }

  double get_square_lower_endpoint(int, int) const override {
    return 0.;
  }

  bool is_nonincreasing_after_lower_endpoint() const override {
    return Varphi::is_nonincreasing_after_lower_endpoint;
  }
  bool is_two_valued() const override {
    return Varphi::is_two_valued;
  }
};

struct Bump_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = true;
  static constexpr bool is_two_valued = false;

  static double set(double radius) {
    return -radius * std::log(2);
  }

  static double apply(double square_distance, double constant) {
    return 1.0 - std::exp(constant / std::sqrt(square_distance));
  }
};

struct Square_bump_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = true;
  static constexpr bool is_two_valued = false;

  static double set(double radius) {
    return -radius * radius * std::log(2);
  }

  static double apply(double square_distance, double constant) {
    return 1.0 - std::exp(constant / square_distance);
  }
};

struct Square_exponential_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = true;
  static constexpr bool is_two_valued = false;

  static double set(double radius) {
    return -std::log(2) / (radius * radius);
  }

  static double apply(double square_distance, double constant) {
    return std::exp(constant * square_distance);
  }
};

struct Exponential_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = true;
  static constexpr bool is_two_valued = false;

  static double set(double radius) {
    return -std::log(2) / radius;
  }

  static double apply(double square_distance, double constant) {
    return std::exp(constant * std::sqrt(square_distance));
  }
};

struct Linear_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = true;
  static constexpr bool is_two_valued = false;

  static double set(double radius) {
    return -1. / radius;
  }

  static double apply(double square_distance, double constant) {
    return std::max<double>(0., 1. + constant * std::sqrt(square_distance));
  }
};

struct Strauss_implementation {
  static constexpr bool is_nonincreasing_after_lower_endpoint = true;
  static constexpr bool is_two_valued = true;

  static double set(double radius) {
    return radius * radius;
  }

  static double apply(double square_distance, double constant) {
    if(square_distance <= constant) {
      return 1.;
    } else {
      return 0.;
    }
  }
};

using Exponential = Short_range_potential<Exponential_implementation>;
using Square_exponential = Short_range_potential<Square_exponential_implementation>;
using Bump = Short_range_potential<Bump_implementation>;
using Square_bump = Short_range_potential<Square_bump_implementation>;
using Strauss = Short_range_potential<Strauss_implementation>;
using Linear = Short_range_potential<Linear_implementation>;

} // namespace potentials

const constexpr char* const short_range_models[] = {
  "exponential",
  "square_exponential",
  "bump",
  "square_bump",
  "Geyer",
  "linear"
};

inline std::shared_ptr<const Potential> make_short_range_object(Rcpp::CharacterVector model,
                                                                Rcpp::NumericMatrix radius) {
  const auto model_string(model[0]);
  if(model_string == short_range_models[0]) {
    return std::make_shared<potentials::Exponential>(radius);
  } else if(model_string == short_range_models[1]) {
    return std::make_shared<potentials::Square_exponential>(radius);
  } else if(model_string == short_range_models[2]) {
    return std::make_shared<potentials::Bump>(radius);
  } else if(model_string == short_range_models[3]) {
    return std::make_shared<potentials::Square_bump>(radius);
  } else if(model_string == short_range_models[4]) {
    return std::make_shared<potentials::Strauss>(radius);
  } else if(model_string == short_range_models[5]) {
    return std::make_shared<potentials::Linear>(radius);
  } else {
    Rcpp::stop("Incorrect model entered. A call to show_short_range_models() will show you the available choices.\n");
  }
}

} // namespace ppjsdm

#endif // INCLUDE_SHORT_RANGE_POTENTIALS
