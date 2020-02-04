#ifndef INCLUDE_PPJSDM_RIPPP
#define INCLUDE_PPJSDM_RIPPP

#include <Rcpp.h>
#include <Rinternals.h>

#include "rppp_single.hpp"

namespace ppjsdm {

template<typename Configuration, typename Window, typename Lambda, typename UpperBound>
inline auto rippp(const Window& window, const Lambda& normalized_lambda, const UpperBound& upper_bound, R_xlen_t number_types) {
  const auto volume(window.volume());

  std::vector<R_xlen_t> number_points(number_types);
  for(R_xlen_t j(0); j < number_types; ++j) {
    const auto points_to_add(R::rpois(volume * static_cast<double>(upper_bound[j])));
    number_points[j] = points_to_add;
  }

  Configuration configuration{};
  // TODO: Might want to .reserve() if Configuration has the member function.
  using size_t = size_t<Configuration>;

  for(R_xlen_t j(0); j < number_types; ++j) {
    const auto points_to_add = static_cast<size_t>(number_points[j]);
    for(size_t i(0); i < points_to_add; ++i) {
      const auto sample(window.sample(j));
      if(normalized_lambda(sample) < 0 || normalized_lambda(sample) > 1) {
        Rcpp::stop("Did not correctly normalise lambda");
      }
      // TODO: As usual, might be able to take log on both sides...
      if(unif_rand() < normalized_lambda(sample)) {
        configuration.push_back(std::move(sample));
      }
    }
  }
  return configuration;
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_RIPPP
