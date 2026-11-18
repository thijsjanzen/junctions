//  Copyright (c) 2026, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

// bare minimal tbb-stub header.
// just enough to to let *this* package pass the
// RCPP_PARALLEL_USE_TBB shenanigan on Alpine Linux

#pragma once

#include <cstdlib>
#include <utility>
#include <Rcpp.h>
#include <RcppParallel.h>   // pull RCPP_PARALLEL_USE_TBB


#if RCPP_PARALLEL_USE_TBB == 0

// everything looks so single-threaded here :(

#include <algorithm>

namespace tbb {

template<typename InputIterator, typename Body>
inline void parallel_for_each(InputIterator first, InputIterator last,
                              Body&& body) {
  std::for_each(first, last, std::forward<Body>(body));
}


template<typename Index, typename Func>
inline void parallel_for(Index first, Index last, const Func& f) {
  for (; first < last; ++first) {
    f(first);
  }
}


template<typename Index, typename Func>
inline void parallel_for(Index first, Index last, Index step,
                         const Func& f) {
  for (; first < last; first += step) {
    f(first);
  }
}


class task_arena {
 public:
  static const int automatic = -1;
  explicit task_arena(int /* max_concurrency */ = automatic) {}
  template<typename Func>
  void execute(const Func& f) { f(); }
};

}  // namespace tbb


// function name is lying.
inline int get_rcpp_num_threads() {
  return 1;
}


#else  // if RCPP_PARALLEL_USE_TBB == 0


// probably the cleanest way to retrieve RcppParallel's concurrency setting
// set by RcppParallel::setThreadOptions(numThreads).
// setThreadOptions("auto") and setThreadOptions(-1) unset the variable;
// task_arena::automatic (-1) then lets TBB use all available cores.
inline int get_rcpp_num_threads() {
  auto* nt_env = std::getenv("RCPP_PARALLEL_NUM_THREADS");
  int n = (nullptr == nt_env) ? 0 : std::atoi(nt_env);
  return n > 0 ? n : static_cast<int>(tbb::task_arena::automatic);
}


#endif
