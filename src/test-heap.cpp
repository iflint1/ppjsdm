#include <Rcpp.h>
#include <testthat.h>

#include "utility/heap.hpp"

#include <algorithm> // std::make_heap, std::sort
#include <cmath> // std::floor, std::fabs
#include <vector> // std::vector

void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

context("Heap") {
  test_that("Get nth element of a heap") {
    set_seed(1);
    const int max_tested_size(20);
    const int number_replications(10);
    for(int i(0); i < max_tested_size; ++i) {
      for(int replication(0); replication < number_replications; ++replication) {
        std::vector<double> vector(i);
        for(int j(0); j < i; ++j) {
          vector[j] = unif_rand();
        }
        auto sorted_vector(vector);

        // Max-heap
        std::make_heap(vector.begin(), vector.end(), std::less<double>{});
        std::sort(sorted_vector.begin(), sorted_vector.end(), std::greater<double>{});

        if(i >= 1) {
          expect_true(ppjsdm::get_nth<0>(vector, std::less<double>{}) == sorted_vector[0]);
        }
        if(i >= 2) {
          expect_true(ppjsdm::get_nth<1>(vector, std::less<double>{}) == sorted_vector[1]);
        }
        if(i >= 3) {
          expect_true(ppjsdm::get_nth<2>(vector, std::less<double>{}) == sorted_vector[2]);
        }

        // Min-heap
        std::make_heap(vector.begin(), vector.end(), std::greater<double>{});
        std::sort(sorted_vector.begin(), sorted_vector.end(), std::less<double>{});

        if(i >= 1) {
          expect_true(ppjsdm::get_nth<0>(vector, std::greater<double>{}) == sorted_vector[0]);
        }
        if(i >= 2) {
          expect_true(ppjsdm::get_nth<1>(vector, std::greater<double>{}) == sorted_vector[1]);
        }
        if(i >= 3) {
          expect_true(ppjsdm::get_nth<2>(vector, std::greater<double>{}) == sorted_vector[2]);
        }
      }
    }
  }
}
