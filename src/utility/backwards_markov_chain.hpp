#ifndef INCLUDE_PPJSDM_BACKWARDS_MARKOV_CHAIN
#define INCLUDE_PPJSDM_BACKWARDS_MARKOV_CHAIN

#include <Rinternals.h>

#include "../configuration/configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"
#include "../simulation/inhomogeneous_ppp.hpp"

#include <tuple> // std::tuple
#include <vector> // std::vector

namespace ppjsdm {
namespace detail {

constexpr double do_not_need_mark = -1.0;

} // namespace detail

template <typename Configuration, typename Model>
class Backwards_Markov_chain {
public:
  explicit Backwards_Markov_chain(const Model& model, R_xlen_t number_types) :
    model_(model),
    last_configuration_(simulate_inhomogeneous_ppp<Configuration>(model, number_types)),
    number_types_(number_types),
    intensity_integral_(model.get_integral()),
    chain_{} {

    }

  auto size() const {
    return chain_.size();
  }

  auto extend_until_T0() {
    using size_t = decltype(ppjsdm::size(last_configuration_));
    const auto initial_last_size(ppjsdm::size(last_configuration_));
    if(initial_last_size == 0) {
      return static_cast<typename decltype(chain_)::size_type>(0);
    }
    Configuration points_not_in_last{};
    while(true) {
      // Do the computations by blocks of size `initial_last_size`, reserving space each time
      chain_.reserve(initial_last_size);
      for(size_t i(0); i < initial_last_size; ++i) {
        const double sum_sizes(ppjsdm::size(points_not_in_last) + ppjsdm::size(last_configuration_));
        const auto v(unif_rand() * (intensity_integral_ + sum_sizes) - intensity_integral_); // Uniformly distributed on [-beta, s_1 + s_2].
        if(v < 0.0) { // Happens with probability beta / (beta + s_1 + s_2).
          insert_uniform_point_in_configuration_and_update_chain(points_not_in_last);
        } else if(v < static_cast<double>(ppjsdm::size(last_configuration_))) { // Happens with probability s_2 / (s_1 + s_2)
          delete_random_point_in_configuration_and_update_chain(last_configuration_);
          if(empty(last_configuration_)) {
            last_configuration_ = points_not_in_last;
            return chain_.size();
          }
        } else {
          delete_random_point_in_configuration_and_update_chain(points_not_in_last);
        }
      }
      R_CheckUserInterrupt();
    }
  }

  template<typename IntegerType>
  void extend_backwards(IntegerType number_extensions) {
    chain_.reserve(number_extensions);
    for(IntegerType i(0); i < number_extensions; ++i) {
      if(unif_rand() * intensity_integral_ + static_cast<double>(ppjsdm::size(last_configuration_)) < intensity_integral_) {
        insert_uniform_point_in_configuration_and_update_chain(last_configuration_);
      } else {
        delete_random_point_in_configuration_and_update_chain(last_configuration_);
      }
    }
  }

  auto get_last_configuration() const {
    return last_configuration_;
  }

  // This is an efficient implementation of the generic CFTP algorithm from Moller and Waagepetersen's book, cf. p. 230 therein.
  // The point process is neither assumed to be attractive nor repulsive, see also p. 361 of ``A primer on perfect simulation
  // for spatial point processes'' by Berthelsen and Moller, and ``Perfect simulation of spatial
  // point processes using dominated coupling from the past with application to a multiscale area-interaction point process''
  // by Ambler and Silverman for a similar (but less general) idea.
  inline auto compute_LU_and_check_coalescence() const {
    auto points_not_in_L(last_configuration_); // U starts with the end value of the chain.
    Configuration L{}; // L is an empty point process.
    const auto chain_size(chain_.size());
    if(chain_size != 0) {
      for(auto n(static_cast<long long int>(chain_size) - 1); n >= 0; --n) {
        const auto& current(chain_[static_cast<std::size_t>(n)]);
        const auto& exp_mark(std::get<0>(current));
        const auto& point(std::get<1>(current));
        if(exp_mark >= 0.0) { // birth
          const auto log_alpha_max(model_.compute_log_alpha_max(point, L, points_not_in_L));
          if(log_alpha_max > 0) {
            Rcpp::stop("Did not correctly normalize the Papangelou intensity");
          }
          if(log_alpha_max + exp_mark > 0) {
            const auto log_alpha_min(model_.compute_log_alpha_min(point, L, points_not_in_L));
            if(log_alpha_min > 0) {
              Rcpp::stop("Did not correctly normalize the Papangelou intensity");
            }
            add_point(log_alpha_min + exp_mark > 0 ? L : points_not_in_L, std::forward<decltype(point)>(point));
          }
        } else { // death
          if(!remove_point(points_not_in_L, point)) {
            remove_point(L, std::forward<decltype(point)>(point));
          }
        }
      }
    }
    return std::pair<bool, Configuration>(empty(points_not_in_L), L);
  }

private:
  Model model_;
  Configuration last_configuration_;
  R_xlen_t number_types_;
  double intensity_integral_;
  std::vector<std::tuple<double, Marked_point>> chain_;

  void insert_uniform_point_in_configuration_and_update_chain(Configuration& configuration) {
    const auto random_type(Rcpp::sample(number_types_, 1, false, R_NilValue, false)[0]);
    const auto point(model_.sample_point(random_type));
    chain_.emplace_back(detail::do_not_need_mark, point);
    add_point(configuration, std::move(point));
  }

  void delete_random_point_in_configuration_and_update_chain(Configuration& configuration) {
    const auto point(remove_random_point(configuration));
    chain_.emplace_back(exp_rand(), std::move(point));
  }
};

template<typename Configuration, typename Model, typename... Args>
inline auto make_backwards_markov_chain(Model&& model, Args&&... args) {
  return Backwards_Markov_chain<Configuration, std::remove_cv_t<std::remove_reference_t<decltype(model)>>>(std::forward<Model>(model), std::forward<Args>(args)...);
}

}  // namespace ppjsdm

#endif  // INCLUDE_PPJSDM_BACKWARDS_MARKOV_CHAIN
