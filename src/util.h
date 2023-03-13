#pragma once

#include "Fish.h"
#include "random_functions.h"

#include <RcppParallel.h>
#include <Rcpp.h>
using namespace Rcpp;


template <typename FISH>
void update_pop(const std::vector<FISH>& old_pop,
                std::vector<FISH>& pop,
                int popSize,
                const std::vector<double>& numRecombinations,
                int num_threads) {

  if (num_threads == 1) {
    rnd_t rndgen;
    for (unsigned i = 0; i < popSize; ++i) {
      int index1 = rndgen.random_number(popSize);
      int index2 = rndgen.random_number(popSize);
      while(index2 == index1) index2 = rndgen.random_number(popSize);

      pop[i] = mate(old_pop[index1], old_pop[index2], numRecombinations,
                        rndgen);
    }
  } else {

    tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);
    rnd_t rndgen;

    tbb::parallel_for(
      tbb::blocked_range<unsigned>(0, popSize),
      [&](const tbb::blocked_range<unsigned>& r) {

        thread_local int seed = get_seed();
        thread_local rnd_t rndgen2(seed);

        for (unsigned i = r.begin(); i < r.end(); ++i) {
          int index1 = rndgen2.random_number(popSize);
          int index2 = rndgen2.random_number(popSize);
          while(index2 == index1) index2 = rndgen2.random_number(popSize);

          pop[i] = mate(old_pop[index1], old_pop[index2], numRecombinations,
                        rndgen2);
        }
      }
    );
  }
}