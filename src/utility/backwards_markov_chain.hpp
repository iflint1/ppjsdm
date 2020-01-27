#ifndef INCLUDE_PPJSDM_BACKWARDS_MARKOV_CHAIN
#define INCLUDE_PPJSDM_BACKWARDS_MARKOV_CHAIN

#include "../configuration/configuration_manipulation.hpp"

namespace ppjsdm {

template <typename Configuration>
class Backwards_Markov_chain {
public:
  explicit Backwards_Markov_chain(Configuration pp) : last_configuration_(std::move(pp)) {}

  auto size() const {
    return (chain_.size() + 1);
  }

  template<typename Window, typename U>
  std::size_t extend_until_T0(U beta, Window& window, R_xlen_t point_types) {
    if(empty(last_configuration_)) {
      return 0;
    }
    Configuration points_not_in_last{};
    std::size_t T0(1);
    while(true) {
      if (unif_rand() * (beta + static_cast<decltype(beta)>(size(points_not_in_last) + size(last_configuration_))) <= beta) {
        insert_uniform_point_in_configuration_and_update_chain(window, points_not_in_last, point_types);
      } else {
        if(unif_rand() * size(last_configuration_) + size(points_not_in_last) < size(last_configuration_)) {  // random point chosen in last_configuration_
          delete_random_point_in_configuration_and_update_chain(last_configuration_);
          if(empty(last_configuration_)) {
            last_configuration_ = points_not_in_last;
            return T0;
          }
        } else {
          delete_random_point_in_configuration_and_update_chain(points_not_in_last);
        }
      }
      ++T0;
    }
  }

  template <typename Window, typename I, typename U>
  void extend_backwards(I number_extensions, U beta, Window& window, R_xlen_t point_types) {
    while(number_extensions--) {
      if(unif_rand() * (beta + static_cast<decltype(beta)>(size(last_configuration_))) <= beta) {
        insert_uniform_point_in_configuration_and_update_chain(window, last_configuration_, point_types);
      } else {
        delete_random_point_in_configuration_and_update_chain(last_configuration_);
      }
    }
  }

  auto get_last_configuration() const {
    return last_configuration_;
  }

  template <typename Function>
  void iterate_forward_in_time(const Function& f) {
    for (auto n(static_cast<long>(chain_.size()) - 1); n >= 0; --n) {
      const auto current(chain_[static_cast<std::size_t>(n)]);
      f(std::get<1>(current), std::get<0>(current), std::get<2>(current));
    }
  }

private:
  Configuration last_configuration_;
  std::vector<std::tuple<bool, P, T>> chain_;

  template <typename Window, typename Configuration>
  void insert_uniform_point_in_configuration_and_update_chain(const Window& window, Configuration& configuration, R_xlen_t point_types) {
    const auto random_type(Rcpp::sample(point_types, 1, false, R_NilValue, false)[0]);
    const auto point(window.sample(random_type));
    chain_.emplace_back(false, point, 0);
    add_point(configuration, std::move(point));
  }

  template <typename Generator, typename Configuration>
  void delete_random_point_in_configuration_and_update_chain(Configuration& configuration) {
    const auto point(remove_random_point(configuration));
    chain_.emplace_back(true, std::move(point), unif_rand());
  }
};

}  // namespace ppjsdm

#endif  // INCLUDE_PPJSDM_BACKWARDS_MARKOV_CHAIN
