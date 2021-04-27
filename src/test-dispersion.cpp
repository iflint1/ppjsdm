#include <Rcpp.h>
#include <testthat.h>

#include "configuration/configuration_manipulation.hpp"
#include "point/point_manipulation.hpp"
#include "point/point_manipulation.hpp"
#include "saturated_model/compute_dispersion.hpp"
#include "saturated_model/compute_dispersion_fitting.hpp"
#include "saturated_model/compute_dispersion_vcov.hpp"
#include "saturated_model/saturated_model.hpp"
#include "saturated_model/potentials/medium_range_potentials.hpp"
#include "saturated_model/potentials/short_range_potentials.hpp"
#include "utility/flatten_strict_upper_triangular.hpp"

#include <cmath> // std::sin
#include <vector> // std::vector

namespace detail {

// Quick and dirty PRNG
long discrete(long& l) {
  const long m = 4294967296;
  const long a = 1103515245;
  const long c = 12345;
  l = (a * l + c) % m;
  return l;
}

double runif(long& l) {
  const long m = 4294967296;
  const long a = 1103515245;
  const long c = 12345;
  l = (a * l + c) % m;
  return static_cast<double>(l) / static_cast<double>(m);
}

} // namespace detail

context("Dispersion") {
  test_that("Compute dispersion by different methods") {
    const int min_tested_size(0);
    const int max_tested_size(11);

    const int min_saturation(2);
    const int max_saturation(3);

    const int ntypes(3);
    const int nmedium(9);

    long state(0);

    using Configuration = std::vector<ppjsdm::Marked_point>;
    for(int i(min_tested_size); i <= max_tested_size; ++i) {
      for(int saturation(min_saturation); saturation <= max_saturation; ++saturation) {
        for(int medium_range_index(0); medium_range_index < nmedium; ++medium_range_index) {
          Configuration configuration(i);
          Configuration other_configuration(i);
          for(int j(0); j < i; ++j) {
            configuration[j] = ppjsdm::Marked_point(detail::runif(state),
                                                    detail::runif(state),
                                                    detail::discrete(state) % ntypes,
                                                    1.0);
            other_configuration[j] = ppjsdm::Marked_point(detail::runif(state),
                                                          detail::runif(state),
                                                          detail::discrete(state) % ntypes,
                                                          1.0);
          }
          Rcpp::NumericMatrix medium_range(ntypes, ntypes);
          Rcpp::NumericMatrix long_range(ntypes, ntypes);
          for(int k1(0); k1 < ntypes; ++k1) {
            for(int k2(k1); k2 < ntypes; ++k2) {
              medium_range(k1, k2) = 0.2 * detail::runif(state);
              medium_range(k2, k1) = medium_range(k1, k2);

              long_range(k1, k2) = medium_range(k1, k2) + 0.2 * detail::runif(state);
              long_range(k2, k1) = long_range(k1, k2);
            }
          }
          const ppjsdm::Saturated_model model(ppjsdm::medium_range_models[medium_range_index], medium_range, long_range, saturation);
          const auto fitting_dispersion(ppjsdm::compute_dispersion_for_fitting(model,
                                                                               ntypes,
                                                                               configuration,
                                                                               other_configuration));
          const auto vcov_dispersion(ppjsdm::compute_dispersion_for_vcov(model,
                                                                         ntypes,
                                                                         configuration));

          for(int type(0); type < ntypes; ++type) {
            for(Configuration::size_type j(0); j < configuration.size(); ++j) {
              const auto dispersion1(ppjsdm::compute_dispersion(model,
                                                                configuration[j],
                                                                             ntypes,
                                                                             configuration)[type]);
              expect_true(exp(dispersion1) == Approx(exp(fitting_dispersion[j][type])));
            }
            for(Configuration::size_type j(0); j < other_configuration.size(); ++j) {
              const auto dispersion1(ppjsdm::compute_dispersion(model,
                                                                other_configuration[j],
                                                                                   ntypes,
                                                                                   configuration)[type]);
              expect_true(exp(dispersion1) == Approx(exp(fitting_dispersion[j + configuration.size()][type])));
            }

            for(Configuration::size_type k1(0); k1 < configuration.size(); ++k1) {
              for(Configuration::size_type k2(k1 + 1); k2 < configuration.size(); ++k2) {
                Configuration configuration_minus(configuration);
                ppjsdm::remove_point(configuration_minus, configuration[k1]);
                ppjsdm::remove_point(configuration_minus, configuration[k2]);
                const auto dispersion1(ppjsdm::compute_dispersion(model,
                                                                  configuration[k1],
                                                                               ntypes,
                                                                               configuration_minus)[type]);
                expect_true(exp(dispersion1) == Approx(exp(vcov_dispersion.first[ppjsdm::encode_linear(k1, k2, configuration.size())][type])));

                const auto dispersion2(ppjsdm::compute_dispersion(model,
                                                                  configuration[k2],
                                                                               ntypes,
                                                                               configuration_minus)[type]);
                expect_true(exp(dispersion2) == Approx(exp(vcov_dispersion.second[ppjsdm::encode_linear(k1, k2, configuration.size())][type])));
              }
            }
          }
        }
      }
    }
  }
}
