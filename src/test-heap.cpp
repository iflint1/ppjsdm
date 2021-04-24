#include <Rcpp.h>
#include <testthat.h>

#include "utility/heap.hpp"

#include <algorithm> // std::make_heap, std::sort
#include <cmath> // std::floor, std::fabs
#include <vector> // std::vector

context("Heap") {
  test_that("Get nth element of a heap") {
    const int max_tested_size(20);
    const int number_replications(10);
    for(int i(1); i < max_tested_size; ++i) {
      for(int replication(0); replication < number_replications; ++replication) {
        std::vector<double> vector(i);
        for(int j(0); j < i; ++j) {
          vector[j] = std::sin(j + replication);
        }
        auto sorted_vector(vector);

        // Max-heap
        std::make_heap(vector.begin(), vector.end(), std::less<double>{});
        std::sort(sorted_vector.begin(), sorted_vector.end(), std::greater<double>{});

        expect_true(ppjsdm::get_nth<0>(vector, std::less<double>{}) == sorted_vector[0]);
        if(i >= 2) {
          expect_true(ppjsdm::get_nth<1>(vector, std::less<double>{}) == sorted_vector[1]);
        }
        if(i >= 3) {
          expect_true(ppjsdm::get_nth<2>(vector, std::less<double>{}) == sorted_vector[2]);
        }

        // Min-heap
        std::make_heap(vector.begin(), vector.end(), std::greater<double>{});
        std::sort(sorted_vector.begin(), sorted_vector.end(), std::less<double>{});

        expect_true(ppjsdm::get_nth<0>(vector, std::greater<double>{}) == sorted_vector[0]);
        if(i >= 2) {
          expect_true(ppjsdm::get_nth<1>(vector, std::greater<double>{}) == sorted_vector[1]);
        }
        if(i >= 3) {
          expect_true(ppjsdm::get_nth<2>(vector, std::greater<double>{}) == sorted_vector[2]);
        }
      }
    }
  }

  test_that("Get nth smallest element of a Max-heap") {
    const int max_tested_size(20);
    const int number_replications(20);
    for(int i(1); i < max_tested_size; ++i) {
      for(int replication(0); replication < number_replications; ++replication) {
        std::vector<double> vector(i);
        for(int j(0); j < i; ++j) {
          vector[j] = std::sin(j + replication);
        }
        auto sorted_vector(vector);

        std::make_heap(vector.begin(), vector.end(), std::less<double>{});
        std::sort(sorted_vector.begin(), sorted_vector.end(), std::greater<double>{});

        // The size of the heap is i, and so the i-th smallest is always the largest.
        // This observation is true whenever the first argument is i.
        expect_true(ppjsdm::get_nth_smallest_if<0>(i - 0, vector, std::less<double>{}, [](const auto){ return false; }) == sorted_vector[0]);
        expect_true(ppjsdm::get_nth_smallest_if<1>(i - 0, vector, std::less<double>{}, [](const auto){ return false; }) == sorted_vector[0]);
        expect_true(ppjsdm::get_nth_smallest_if<2>(i - 0, vector, std::less<double>{}, [](const auto){ return false; }) == sorted_vector[0]);
        if(i >= 2) {
          // The size of the heap is i, and we want the 'i-1'-th smallest, i.e. the second-largest.
          expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, std::less<double>{}, [](const auto){ return false; }) == sorted_vector[1]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, std::less<double>{}, [](const auto){ return false; }) == sorted_vector[1]);

          // The cases below correspond to the removal of sorted_vector[l].
          // If we remove the largest one, the 'i-1'-th smallest does not change.
          // In other cases, it becomes the largest one.
          expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, std::less<double>{}, [&sorted_vector](const auto element){ return element >= sorted_vector[0]; }) == sorted_vector[1]);
          expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, std::less<double>{}, [&sorted_vector](const auto element){ return element >= sorted_vector[1]; }) == sorted_vector[0]);
          if(i >= 3) {
            expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, std::less<double>{}, [&sorted_vector](const auto element){ return element >= sorted_vector[2]; }) == sorted_vector[0]);
          }
          if(i >= 4) {
            expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, std::less<double>{}, [&sorted_vector](const auto element){ return element >= sorted_vector[3]; }) == sorted_vector[0]);
          }
        }
        if(i >= 3) {
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, std::less<double>{}, [](const auto){ return false; }) == sorted_vector[2]);

          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, std::less<double>{}, [&sorted_vector](const auto element){ return element >= sorted_vector[0]; }) == sorted_vector[1]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, std::less<double>{}, [&sorted_vector](const auto element){ return element >= sorted_vector[1]; }) == sorted_vector[0]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, std::less<double>{}, [&sorted_vector](const auto element){ return element >= sorted_vector[2]; }) == sorted_vector[0]);

          // The cases below correspond to the removal of sorted_vector[l].
          // If we remove the largest one or second-largest one, the 'i-1'-th smallest does not change.
          // In other cases, it becomes the second-largest one instead of the third-largest.
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, std::less<double>{}, [&sorted_vector](const auto element){ return element >= sorted_vector[0]; }) == sorted_vector[2]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, std::less<double>{}, [&sorted_vector](const auto element){ return element >= sorted_vector[1]; }) == sorted_vector[2]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, std::less<double>{}, [&sorted_vector](const auto element){ return element >= sorted_vector[2]; }) == sorted_vector[1]);
          if(i >= 4) {
            expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, std::less<double>{}, [&sorted_vector](const auto element){ return element >= sorted_vector[3]; }) == sorted_vector[0]);
            expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, std::less<double>{}, [&sorted_vector](const auto element){ return element >= sorted_vector[3]; }) == sorted_vector[1]);
          }
        }
      }
    }
  }

  test_that("Get nth largest element of a Min-heap") {
    const int max_tested_size(20);
    const int number_replications(20);
    for(int i(1); i < max_tested_size; ++i) {
      for(int replication(0); replication < number_replications; ++replication) {
        std::vector<double> vector(i);
        for(int j(0); j < i; ++j) {
          vector[j] = std::sin(j + replication);
        }
        auto sorted_vector(vector);

        std::make_heap(vector.begin(), vector.end(), std::greater<double>{});
        std::sort(sorted_vector.begin(), sorted_vector.end(), std::less<double>{});

        expect_true(ppjsdm::get_nth_smallest_if<0>(i - 0, vector, std::greater<double>{}, [](const auto){ return false; }) == sorted_vector[0]);
        expect_true(ppjsdm::get_nth_smallest_if<1>(i - 0, vector, std::greater<double>{}, [](const auto){ return false; }) == sorted_vector[0]);
        expect_true(ppjsdm::get_nth_smallest_if<2>(i - 0, vector, std::greater<double>{}, [](const auto){ return false; }) == sorted_vector[0]);
        if(i >= 2) {
          expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, std::greater<double>{}, [](const auto){ return false; }) == sorted_vector[1]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, std::greater<double>{}, [](const auto){ return false; }) == sorted_vector[1]);

          expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, std::greater<double>{}, [&sorted_vector](const auto element){ return element <= sorted_vector[0]; }) == sorted_vector[1]);
          expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, std::greater<double>{}, [&sorted_vector](const auto element){ return element <= sorted_vector[1]; }) == sorted_vector[0]);
          if(i >= 3) {
            expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, std::greater<double>{}, [&sorted_vector](const auto element){ return element <= sorted_vector[2]; }) == sorted_vector[0]);
          }
          if(i >= 4) {
            expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, std::greater<double>{}, [&sorted_vector](const auto element){ return element <= sorted_vector[3]; }) == sorted_vector[0]);
          }
        }
        if(i >= 3) {
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, std::greater<double>{}, [](const auto){ return false; }) == sorted_vector[2]);

          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, std::greater<double>{}, [&sorted_vector](const auto element){ return element <= sorted_vector[0]; }) == sorted_vector[1]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, std::greater<double>{}, [&sorted_vector](const auto element){ return element <= sorted_vector[1]; }) == sorted_vector[0]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, std::greater<double>{}, [&sorted_vector](const auto element){ return element <= sorted_vector[2]; }) == sorted_vector[0]);

          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, std::greater<double>{}, [&sorted_vector](const auto element){ return element <= sorted_vector[0]; }) == sorted_vector[2]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, std::greater<double>{}, [&sorted_vector](const auto element){ return element <= sorted_vector[1]; }) == sorted_vector[2]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, std::greater<double>{}, [&sorted_vector](const auto element){ return element <= sorted_vector[2]; }) == sorted_vector[1]);
          if(i >= 4) {
            expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, std::greater<double>{}, [&sorted_vector](const auto element){ return element <= sorted_vector[3]; }) == sorted_vector[0]);
            expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, std::greater<double>{}, [&sorted_vector](const auto element){ return element <= sorted_vector[3]; }) == sorted_vector[1]);
          }
        }
      }
    }
  }
}
