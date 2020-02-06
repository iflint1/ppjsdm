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

constexpr double do_not_need_mark = 2.0;

} // namespace detail

template <typename Configuration, typename Intensity>
class Backwards_Markov_chain {
public:
  explicit Backwards_Markov_chain(const Intensity& intensity, R_xlen_t number_types) :
    intensity_(intensity),
    last_configuration_(simulate_inhomogeneous_ppp<Configuration>(intensity, number_types)) {}

  template<typename Window>
  auto extend_until_T0(double beta, const Window& window, R_xlen_t number_types) {
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
        const auto beta_plus_sizes(beta + static_cast<double>(ppjsdm::size(points_not_in_last) + ppjsdm::size(last_configuration_)));
        const auto v(unif_rand() * beta_plus_sizes - beta); // Uniformly distributed on [-beta, s_1 + s_2].
        if(v < 0.0) { // Happens with probability beta / (beta + s_1 + s_2).
          insert_uniform_point_in_configuration_and_update_chain(window, points_not_in_last, number_types);
        } else {
          const auto w(v / static_cast<double>(ppjsdm::size(last_configuration_))); // Uniformly distributed on [0, 1 + s_1 / s_2].
          if(w < 1.0) { // Happens with probability s_2 / (s_1 + s_2)
            delete_random_point_in_configuration_and_update_chain(last_configuration_, w);
            if(empty(last_configuration_)) {
              last_configuration_ = points_not_in_last;
              return chain_.size();
            }
          } else {
            const auto x((w - 1.0) * static_cast<double>(ppjsdm::size(last_configuration_)) / static_cast<double>(ppjsdm::size(points_not_in_last)));
            delete_random_point_in_configuration_and_update_chain(points_not_in_last, x);
          }
        }
      }
      R_CheckUserInterrupt();
    }
  }

  template<typename Window, typename IntegerType>
  void extend_backwards(IntegerType number_extensions, double beta, const Window& window, R_xlen_t number_types) {
    chain_.reserve(number_extensions);
    for(IntegerType i(0); i < number_extensions; ++i) {
      const auto beta_plus_size(beta + static_cast<double>(ppjsdm::size(last_configuration_)));
      const auto v(unif_rand() * beta_plus_size - beta); // Uniformly distributed on [-beta, s].
      if(v < 0.0) { // Happens with probability beta / (beta + s).
        insert_uniform_point_in_configuration_and_update_chain(window, last_configuration_, number_types);
      } else {
        const auto w(v / static_cast<double>(ppjsdm::size(last_configuration_))); // Uniformly distributed on [0, 1].
        delete_random_point_in_configuration_and_update_chain(last_configuration_, w);
      }
    }
  }

  auto get_last_configuration() const {
    return last_configuration_;
  }

  template<typename F, typename G>
  void iterate_forward_in_time(const F& birth, const G& death) {
    const auto chain_size(chain_.size());
    if(chain_size == 0) {
      return;
    } else {
      for(auto n(static_cast<long long int>(chain_size) - 1); n >= 0; --n) {
        const auto& current(chain_[static_cast<std::size_t>(n)]);
        if(std::get<0>(current) <= 1.0) {
          birth(std::get<1>(current), std::get<0>(current));
        } else {
          death(std::get<1>(current));
        }
      }
    }
  }

private:
  Intensity intensity_;
  Configuration last_configuration_;
  std::vector<std::tuple<double, Marked_point>> chain_;

  template<typename Window>
  void insert_uniform_point_in_configuration_and_update_chain(const Window& window, Configuration& configuration, R_xlen_t number_types) {
    const auto random_type(Rcpp::sample(number_types, 1, false, R_NilValue, false)[0]);
    const auto point(window.sample(random_type));
    chain_.emplace_back(detail::do_not_need_mark, point);
    add_point(configuration, std::move(point));
  }

  void delete_random_point_in_configuration_and_update_chain(Configuration& configuration, double uniform_mark) {
    const auto point(remove_random_point(configuration));
    chain_.emplace_back(uniform_mark, std::move(point));
  }
};

template<typename Configuration, typename Intensity, typename... Args>
inline auto make_backwards_markov_chain(Intensity&& intensity, Args&&... args) {
  return Backwards_Markov_chain<Configuration, std::remove_cv_t<std::remove_reference_t<decltype(intensity)>>>(std::forward<Intensity>(intensity), std::forward<Args>(args)...);
}

}  // namespace ppjsdm

#endif  // INCLUDE_PPJSDM_BACKWARDS_MARKOV_CHAIN
