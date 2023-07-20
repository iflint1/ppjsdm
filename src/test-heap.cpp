#include <Rcpp.h>
#include <testthat.h>

#include "utility/heap.hpp"

#include <algorithm> // std::make_heap, std::sort
#include <cmath> // std::sin
#include <vector> // std::vector

namespace detail {

struct Acc {
  template<typename S, typename T>
  auto operator()(S count, T element) const { return count + element; }
};

} // namespace detail

context("Heap") {
  test_that("Get nth element of a heap") {
    const int max_tested_size(10);
    const int number_replications(5);
    for(int i(0); i < max_tested_size; ++i) {
      for(int replication(0); replication < number_replications; ++replication) {
        std::vector<double> vector(i);
        for(int j(0); j < i; ++j) {
          vector[j] = std::sin(j + replication);
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

  test_that("Get nth smallest element of a Max-heap") {
    const int max_tested_size(10);
    const int number_replications(5);
    for(int i(1); i < max_tested_size; ++i) {
      for(int replication(0); replication < number_replications; ++replication) {
        std::vector<double> vector(i);
        for(int j(0); j < i; ++j) {
          vector[j] = std::sin(j + replication);
        }
        auto sorted_vector(vector);

        using Greater = std::greater<double>;
        using Less = std::less<double>;

        std::make_heap(vector.begin(), vector.end(), Less{});
        std::sort(sorted_vector.begin(), sorted_vector.end(), Greater{});

        // The size of the heap is i, and so the i-th smallest is always the largest.
        // This observation is true whenever the first argument is i.
        expect_true(ppjsdm::get_nth_smallest_if<0>(i - 0, vector, [](const auto){ return false; }, Less{}) == sorted_vector[0]);
        expect_true(ppjsdm::get_nth_smallest_if<1>(i - 0, vector, [](const auto){ return false; }, Less{}) == sorted_vector[0]);
        expect_true(ppjsdm::get_nth_smallest_if<2>(i - 0, vector, [](const auto){ return false; }, Less{}) == sorted_vector[0]);
        if(i >= 2) {
          // The size of the heap is i, and we want the 'i-1'-th smallest, i.e. the second-largest.
          expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, [](const auto){ return false; }, Less{}) == sorted_vector[1]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, [](const auto){ return false; }, Less{}) == sorted_vector[1]);

          // The cases below correspond to the removal of sorted_vector[l].
          // If we remove the largest one, the 'i-1'-th smallest does not change.
          // In other cases, it becomes the largest one.
          expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, [&sorted_vector](const auto element){ return element >= sorted_vector[0]; }, Less{}) == sorted_vector[1]);
          expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, [&sorted_vector](const auto element){ return element >= sorted_vector[1]; }, Less{}) == sorted_vector[0]);
          if(i >= 3) {
            expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, [&sorted_vector](const auto element){ return element >= sorted_vector[2]; }, Less{}) == sorted_vector[0]);
          }
          if(i >= 4) {
            expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, [&sorted_vector](const auto element){ return element >= sorted_vector[3]; }, Less{}) == sorted_vector[0]);
          }
        }
        if(i >= 3) {
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, [](const auto){ return false; }, Less{}) == sorted_vector[2]);

          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, [&sorted_vector](const auto element){ return element >= sorted_vector[0]; }, Less{}) == sorted_vector[1]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, [&sorted_vector](const auto element){ return element >= sorted_vector[1]; }, Less{}) == sorted_vector[0]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, [&sorted_vector](const auto element){ return element >= sorted_vector[2]; }, Less{}) == sorted_vector[0]);

          // The cases below correspond to the removal of sorted_vector[l].
          // If we remove the largest one or second-largest one, the 'i-1'-th smallest does not change.
          // In other cases, it becomes the second-largest one instead of the third-largest.
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, [&sorted_vector](const auto element){ return element >= sorted_vector[0]; }, Less{}) == sorted_vector[2]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, [&sorted_vector](const auto element){ return element >= sorted_vector[1]; }, Less{}) == sorted_vector[2]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, [&sorted_vector](const auto element){ return element >= sorted_vector[2]; }, Less{}) == sorted_vector[1]);
          if(i >= 4) {
            expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, [&sorted_vector](const auto element){ return element >= sorted_vector[3]; }, Less{}) == sorted_vector[0]);
            expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, [&sorted_vector](const auto element){ return element >= sorted_vector[3]; }, Less{}) == sorted_vector[1]);
          }
        }
      }
    }
  }

  test_that("Get nth largest element of a Min-heap") {
    const int max_tested_size(10);
    const int number_replications(5);
    for(int i(1); i < max_tested_size; ++i) {
      for(int replication(0); replication < number_replications; ++replication) {
        std::vector<double> vector(i);
        for(int j(0); j < i; ++j) {
          vector[j] = std::sin(j + replication);
        }
        auto sorted_vector(vector);

        using Greater = std::greater<double>;
        using Less = std::less<double>;
        std::make_heap(vector.begin(), vector.end(), Greater{});
        std::sort(sorted_vector.begin(), sorted_vector.end(), Less{});

        expect_true(ppjsdm::get_nth_smallest_if<0>(i - 0, vector, [](const auto){ return false; }, Greater{}) == sorted_vector[0]);
        expect_true(ppjsdm::get_nth_smallest_if<1>(i - 0, vector, [](const auto){ return false; }, Greater{}) == sorted_vector[0]);
        expect_true(ppjsdm::get_nth_smallest_if<2>(i - 0, vector, [](const auto){ return false; }, Greater{}) == sorted_vector[0]);
        if(i >= 2) {
          expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, [](const auto){ return false; }, Greater{}) == sorted_vector[1]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, [](const auto){ return false; }, Greater{}) == sorted_vector[1]);

          expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, [&sorted_vector](const auto element){ return element <= sorted_vector[0]; }, Greater{}) == sorted_vector[1]);
          expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, [&sorted_vector](const auto element){ return element <= sorted_vector[1]; }, Greater{}) == sorted_vector[0]);
          if(i >= 3) {
            expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, [&sorted_vector](const auto element){ return element <= sorted_vector[2]; }, Greater{}) == sorted_vector[0]);
          }
          if(i >= 4) {
            expect_true(ppjsdm::get_nth_smallest_if<1>(i - 1, vector, [&sorted_vector](const auto element){ return element <= sorted_vector[3]; }, Greater{}) == sorted_vector[0]);
          }
        }
        if(i >= 3) {
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, [](const auto){ return false; }, Greater{}) == sorted_vector[2]);

          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, [&sorted_vector](const auto element){ return element <= sorted_vector[0]; }, Greater{}) == sorted_vector[1]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, [&sorted_vector](const auto element){ return element <= sorted_vector[1]; }, Greater{}) == sorted_vector[0]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, [&sorted_vector](const auto element){ return element <= sorted_vector[2]; }, Greater{}) == sorted_vector[0]);

          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, [&sorted_vector](const auto element){ return element <= sorted_vector[0]; }, Greater{}) == sorted_vector[2]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, [&sorted_vector](const auto element){ return element <= sorted_vector[1]; }, Greater{}) == sorted_vector[2]);
          expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, [&sorted_vector](const auto element){ return element <= sorted_vector[2]; }, Greater{}) == sorted_vector[1]);
          if(i >= 4) {
            expect_true(ppjsdm::get_nth_smallest_if<2>(i - 1, vector, [&sorted_vector](const auto element){ return element <= sorted_vector[3]; }, Greater{}) == sorted_vector[0]);
            expect_true(ppjsdm::get_nth_smallest_if<2>(i - 2, vector, [&sorted_vector](const auto element){ return element <= sorted_vector[3]; }, Greater{}) == sorted_vector[1]);
          }
        }
      }
    }
  }

  test_that("Accumulate n smallest points of a Max-heap") {
    const int max_tested_size(10);
    const int number_replications(5);
    for(int i(1); i < max_tested_size; ++i) {
      for(int replication(0); replication < number_replications; ++replication) {
        std::vector<double> vector(i);
        for(int j(0); j < i; ++j) {
          vector[j] = std::sin(j + replication);
        }
        auto sorted_vector(vector);

        using Greater = std::greater<double>;
        using Less = std::less<double>;

        std::make_heap(vector.begin(), vector.end(), Less{});
        std::sort(sorted_vector.begin(), sorted_vector.end(), Greater{});

        expect_true(ppjsdm::accumulate_n_smallest<0>(i, vector, Less{}, detail::Acc{}) == Approx(std::accumulate(sorted_vector.begin(), sorted_vector.end(), 0., detail::Acc{})));
        expect_true(ppjsdm::accumulate_n_smallest<1>(i, vector, Less{}, detail::Acc{}) == Approx(std::accumulate(sorted_vector.begin(), sorted_vector.end(), 0., detail::Acc{})));
        expect_true(ppjsdm::accumulate_n_smallest<2>(i, vector, Less{}, detail::Acc{}) == Approx(std::accumulate(sorted_vector.begin(), sorted_vector.end(), 0., detail::Acc{})));
        if(i > 1) {
          expect_true(ppjsdm::accumulate_n_smallest<1>(i - 1, vector, Less{}, detail::Acc{}) == Approx(std::accumulate(sorted_vector.begin() + 1, sorted_vector.end(), 0., detail::Acc{})));
          expect_true(ppjsdm::accumulate_n_smallest<2>(i - 1, vector, Less{}, detail::Acc{}) == Approx(std::accumulate(sorted_vector.begin() + 1, sorted_vector.end(), 0., detail::Acc{})));
        }
        if(i > 2) {
          expect_true(ppjsdm::accumulate_n_smallest<2>(i - 2, vector, Less{}, detail::Acc{}) == Approx(std::accumulate(sorted_vector.begin() + 2, sorted_vector.end(), 0., detail::Acc{})));
        }
      }
    }
  }

  test_that("Accumulate n Largest points of a Min-heap") {
    const int max_tested_size(10);
    const int number_replications(5);
    for(int i(1); i < max_tested_size; ++i) {
      for(int replication(0); replication < number_replications; ++replication) {
        std::vector<double> vector(i);
        for(int j(0); j < i; ++j) {
          vector[j] = std::sin(j + replication);
        }
        auto sorted_vector(vector);

        using Greater = std::greater<double>;
        using Less = std::less<double>;

        std::make_heap(vector.begin(), vector.end(), Greater{});
        std::sort(sorted_vector.begin(), sorted_vector.end(), Less{});

        expect_true(ppjsdm::accumulate_n_smallest<0>(i, vector, Greater{}, detail::Acc{}) == Approx(std::accumulate(sorted_vector.begin(), sorted_vector.end(), 0., detail::Acc{})));
        expect_true(ppjsdm::accumulate_n_smallest<1>(i, vector, Greater{}, detail::Acc{}) == Approx(std::accumulate(sorted_vector.begin(), sorted_vector.end(), 0., detail::Acc{})));
        expect_true(ppjsdm::accumulate_n_smallest<2>(i, vector, Greater{}, detail::Acc{}) == Approx(std::accumulate(sorted_vector.begin(), sorted_vector.end(), 0., detail::Acc{})));
        if(i > 1) {
          expect_true(ppjsdm::accumulate_n_smallest<1>(i - 1, vector, Greater{}, detail::Acc{}) == Approx(std::accumulate(sorted_vector.begin() + 1, sorted_vector.end(), 0., detail::Acc{})));
          expect_true(ppjsdm::accumulate_n_smallest<2>(i - 1, vector, Greater{}, detail::Acc{}) == Approx(std::accumulate(sorted_vector.begin() + 1, sorted_vector.end(), 0., detail::Acc{})));
        }
        if(i > 2) {
          expect_true(ppjsdm::accumulate_n_smallest<2>(i - 2, vector, Greater{}, detail::Acc{}) == Approx(std::accumulate(sorted_vector.begin() + 2, sorted_vector.end(), 0., detail::Acc{})));
        }
      }
    }
  }
}
