#ifndef INCLUDE_PPJSDM_BACKWARDS_MARKOV_CHAIN
#define INCLUDE_PPJSDM_BACKWARDS_MARKOV_CHAIN

#include <Rinternals.h>

#include "../configuration/configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"
#include "../simulation/inhomogeneous_ppp.hpp"

#include <tuple> // std::tuple
#include <type_traits> // std::remove_cv_t, std::remove_reference_t
#include <utility> // std::forward
#include <vector> // std::vector

namespace ppjsdm {
namespace detail {

constexpr double do_not_need_mark = -1.0;
constexpr double always_in_L = -2.0;

template <typename Configuration, typename Model>
class Backwards_Markov_chain {
public:
  explicit Backwards_Markov_chain(const Model& model) :
    model_(model),
    last_configuration_(simulate_inhomogeneous_ppp<Configuration>(model.get_window(),
                                                                  [&model](const auto& point) { return model.get_log_normalized_bounding_intensity(point); },
                                                                  model.get_upper_bound(),
                                                                  model.get_number_types())),
                                                                  intensity_integral_(model.get_integral()),
                                                                  chain_{} {}

  auto size() const {
    return chain_.size();
  }

  void extend_until_T0() {
    const auto initial_last_size(ppjsdm::size(last_configuration_));
    if(initial_last_size != 0) {
      Configuration points_not_in_last{};
      while(true) {
        // Do the computations by blocks of size `initial_last_size`, reserving space each time
        chain_.reserve(chain_.size() + initial_last_size);
        ppjsdm::reserve_if_possible(points_not_in_last, points_not_in_last.size() + initial_last_size);
        using size_t = decltype(ppjsdm::size(last_configuration_));
        for(size_t i(0); i < initial_last_size; ++i) {
          const double sum_sizes(ppjsdm::size(points_not_in_last) + ppjsdm::size(last_configuration_));
          const auto v(unif_rand() * (intensity_integral_ + sum_sizes) - intensity_integral_); // Uniformly distributed on [-beta, s_1 + s_2].
          if(v < 0.0) { // Happens with probability beta / (beta + s_1 + s_2).
            insert_uniform_point_in_configuration_and_update_chain(points_not_in_last);
          } else if(v < static_cast<double>(ppjsdm::size(last_configuration_))) { // Happens with probability s_2 / (s_1 + s_2).
            delete_random_point_in_configuration_and_update_chain(last_configuration_);
            if(empty(last_configuration_)) {
              last_configuration_ = std::move(points_not_in_last);
              return;
            }
          } else {
            delete_random_point_in_configuration_and_update_chain(points_not_in_last);
          }
        }
        R_CheckUserInterrupt();
      }
    }
  }

  template<typename IntegerType>
  void extend_backwards(IntegerType number_extensions) {
    chain_.reserve(chain_.size() + number_extensions);
    for(IntegerType i(0); i < number_extensions; ++i) {
      if(unif_rand() * (intensity_integral_ + static_cast<double>(ppjsdm::size(last_configuration_))) < intensity_integral_) {
        insert_uniform_point_in_configuration_and_update_chain(last_configuration_);
      } else {
        delete_random_point_in_configuration_and_update_chain(last_configuration_);
      }
    }
  }

  // This is an efficient implementation of the generic CFTP algorithm from Moller and Waagepetersen's book, cf. p. 230 therein.
  // The point process is neither assumed to be attractive nor repulsive, see also p. 361 of ``A primer on perfect simulation
  // for spatial point processes'' by Berthelsen and Moller, and ``Perfect simulation of spatial
  // point processes using dominated coupling from the past with application to a multiscale area-interaction point process''
  // by Ambler and Silverman for a similar (but less general) idea.
  auto compute_LU_and_check_coalescence() const {
    Configuration L_complement(last_configuration_); // U starts with the end value of the chain.
    Configuration L{}; // L is an empty configuration.
    // Reserve a bit extra space for each of the configurations
    ppjsdm::reserve_if_possible(L, 2 * ppjsdm::size(last_configuration_));
    ppjsdm::reserve_if_possible(L_complement, 2 * ppjsdm::size(last_configuration_));
    const auto chain_size(chain_.size());
    if(chain_size != 0) {
      for(auto n(static_cast<long long int>(chain_size) - 1); n >= 0; --n) {
        const auto exp_mark(std::get<0>(chain_[static_cast<std::size_t>(n)]));
        // TODO: Write extensive tests for the functions compute_log_alpha_min_lower_bound, compute_log_alpha_max, etc since I'm not making any checks here.
        if(exp_mark >= 0.0) { // birth
          model_.add_to_L_or_U(exp_mark, std::get<1>(chain_[static_cast<std::size_t>(n)]), L, L_complement);
        } else if(exp_mark == detail::do_not_need_mark) { // death
          const auto& point(std::get<1>(chain_[static_cast<std::size_t>(n)]));
          if(!remove_point(L_complement, point)) {
            remove_point(L, point);
          }
        } else {  // Avoid computations above in this case.
          add_point(L, std::get<1>(chain_[static_cast<std::size_t>(n)]));
        }
      }
    }
    return std::pair<bool, Configuration>(empty(L_complement), L);
  }

private:
  Model model_;
  Configuration last_configuration_;
  double intensity_integral_;
  std::vector<std::tuple<double, Marked_point>> chain_;

  void insert_uniform_point_in_configuration_and_update_chain(Configuration& configuration) {
    const auto point(model_.sample_point_from_bounding_intensity());
    chain_.emplace_back(detail::do_not_need_mark, point);
    add_point(configuration, std::move(point));
  }

  void delete_random_point_in_configuration_and_update_chain(Configuration& configuration) {
    const auto point(remove_random_point(configuration));
    const auto e(exp_rand());
    if(model_.compute_log_alpha_min_lower_bound(get_type(point)) + e > 0) {
      chain_.emplace_back(detail::always_in_L, std::move(point));
    } else {
      chain_.emplace_back(e, std::move(point));
    }
  }
};

} // namespace detail

template<typename Configuration, typename Model, typename... Args>
inline auto make_backwards_markov_chain(Model&& model, Args&&... args) {
  return detail::Backwards_Markov_chain<Configuration, std::remove_cv_t<std::remove_reference_t<decltype(model)>>>(std::forward<Model>(model), std::forward<Args>(args)...);
}

}  // namespace ppjsdm

#endif  // INCLUDE_PPJSDM_BACKWARDS_MARKOV_CHAIN
