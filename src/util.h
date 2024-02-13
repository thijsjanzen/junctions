#pragma once

#include <RcppParallel.h>

inline size_t get_rcpp_num_threads() {
    auto* nt_env = std::getenv("RCPP_PARALLEL_NUM_THREADS");
    return (nullptr == nt_env) 
      ? tbb::task_arena::automatic  // -1
      : static_cast<size_t>(std::atoi(nt_env));
}

inline void set_num_threads() {
    auto num_threads = get_rcpp_num_threads();
    auto global_control = tbb::global_control(tbb::global_control::max_allowed_parallelism, num_threads);
}