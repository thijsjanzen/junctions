#include <stdio.h>
#include <vector>
#include "Output.h"
#include "Fish.h"
#include "random_functions.h"

#include <thread>
#include <chrono>
#include <functional>

#include <RcppParallel.h>
#include <Rcpp.h>
using namespace Rcpp;

void update_pop_multi(const std::vector<Fish_multi>& old_pop,
                std::vector<Fish_multi>& pop,
                int popSize,
                std::vector<double> numRecombinations,
                int num_threads) {

  if (num_threads == 1) {
    rnd_t rndgen;
    for (unsigned i = 0; i < popSize; ++i) {
      int index1 = rndgen.random_number(popSize);
      int index2 = rndgen.random_number(popSize);
      while(index2 == index1) index2 = rndgen.random_number(popSize);

      pop[i] = mate_multi(old_pop[index1], old_pop[index2], numRecombinations,
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

          pop[i] = mate_multi(old_pop[index1], old_pop[index2], numRecombinations,
                              rndgen2);
        }
      }
    );
  }
}

Output simulation_phased_nonphased_multi(int popSize,
                                         double initRatio,
                                         int maxTime,
                                         const std::vector<double>& numRecombinations,
                                         const std::vector< std::vector< double >>&  markers,
                                         const NumericVector& time_points,
                                         bool verbose,
                                         bool record_true_junctions,
                                         int num_indiv_sampled,
                                         int num_threads,
                                         rnd_t& rndgen)    {

  Output O;
  std::vector< Fish_multi > Pop(popSize);

  O.markers_multi = markers;

  Fish_multi parent1 = Fish_multi(0, numRecombinations.size());
  Fish_multi parent2 = Fish_multi(1, numRecombinations.size());

  for (int i = 0; i < popSize; ++i) {
    Fish_multi p1 = parent2;
    Fish_multi p2 = parent2;

    if (rndgen.uniform() < initRatio) {
      p1 = parent1;
    }
    if (rndgen.uniform() < initRatio) {
      p2 = parent1;
    }

    Pop[i] = mate_multi(p1, p2, numRecombinations, rndgen);
  }

  if (verbose) Rcout << "0--------25--------50--------75--------100\n";
  if (verbose) Rcout << "*";
  int updateFreq = maxTime / 20;
  if (updateFreq < 1) updateFreq = 1;

  //  Rcout << "starting simulation\n"; force_output();

  for (size_t t = 0; t <= maxTime; ++t) {
    if (is_in_time_points(t, time_points)) {
      O.update_unphased(Pop, t, record_true_junctions, numRecombinations,
                        num_indiv_sampled);
    }

    std::vector< Fish_multi > newGeneration(popSize);

    update_pop_multi(Pop, newGeneration, popSize, numRecombinations, num_threads);

    Pop.swap(newGeneration);

    if (verbose) {
      if(t % updateFreq == 0) {
        Rcout << "**";
      }


    }
    Rcpp::checkUserInterrupt();
  }
  if (verbose) Rcout << "\n";
  return O;
}

// [[Rcpp::export]]
List sim_phased_unphased_multi_cpp(int pop_size,
                                   double freq_ancestor_1,
                                   int total_runtime,
                                   std::vector<double> size_in_morgan,
                                   NumericMatrix markers,
                                   NumericVector time_points,
                                   bool verbose,
                                   bool record_true_junctions,
                                   int num_indiv_sampled,
                                   int num_threads) {

  rnd_t rndgen;

 // std::vector< double > marker_dist(markers.begin(), markers.end());

  std::vector< std::vector< double >> markers_cpp(markers.nrow(),
                                                  std::vector<double>(markers.ncol()));
  for (size_t i = 0; i < markers.nrow(); ++i) {
    for (size_t j = 0; j < markers.ncol(); ++j) {
      markers_cpp[i][j] = markers(i, j);
    }
  }



  Output O = simulation_phased_nonphased_multi(pop_size,
                                               freq_ancestor_1,
                                               total_runtime,
                                               size_in_morgan,
                                               markers_cpp,
                                               time_points,
                                               verbose,
                                               record_true_junctions,
                                               num_indiv_sampled,
                                               num_threads,
                                               rndgen);

  int num_rows = O.results.size();
  int num_cols = O.results[0].size();

  NumericMatrix output_matrix(num_rows, num_cols);
  for (int i = 0; i < num_rows; ++i) {
    for (int j = 0; j < num_cols; ++j) {
      output_matrix(i,j) = O.results[i][j];
    }
  }

  if (record_true_junctions) {
    int num_rows_t = O.true_results.size();
    int num_cols_t = O.true_results[0].size();
    NumericMatrix output_matrix_true(num_rows_t, num_cols_t);
    for(int i = 0; i < num_rows_t; ++i) {
      for(int j = 0; j < num_cols_t; ++j) {
        output_matrix_true(i,j) = O.true_results[i][j];
      }
    }

    return List::create(Named("results") = output_matrix,
                        Named("true_results") = output_matrix_true);
  }
  return List::create(Named("results") = output_matrix);
}