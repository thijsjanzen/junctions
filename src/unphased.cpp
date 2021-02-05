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

int get_seed();

void update_pop(const std::vector<Fish_inf>& old_pop,
                std::vector<Fish_inf>& pop,
                int popSize,
                int numRecombinations,
                int num_threads) {

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

        pop[i] = mate_inf(old_pop[index1], old_pop[index2], numRecombinations,
                          rndgen2);
      }
    }
  );
}


Output simulation_phased_nonphased(int popSize,
                                   double initRatio,
                                   int maxTime,
                                   double numRecombinations,
                                   std::vector< double >  markers,
                                   const NumericVector& time_points,
                                   bool verbose,
                                   bool record_true_junctions,
                                   int num_indiv_sampled,
                                   int num_threads,
                                   rnd_t& rndgen)    {

  Output O;
  std::vector< Fish_inf > Pop(popSize);

  O.markers = markers;

  Fish_inf parent1 = Fish_inf(0);
  Fish_inf parent2 = Fish_inf(1);

  for (int i = 0; i < popSize; ++i) {
    Fish_inf p1 = parent2;
    Fish_inf p2 = parent2;

    if (rndgen.uniform() < initRatio) {
      p1 = parent1;
    }
    if (rndgen.uniform() < initRatio) {
      p2 = parent1;
    }

    Pop[i] = mate_inf(p1,p2, numRecombinations, rndgen);
  }

  if (verbose) Rcout << "0--------25--------50--------75--------100\n";
  if (verbose) Rcout << "*";
  int updateFreq = maxTime / 20;
  if (updateFreq < 1) updateFreq = 1;

//  Rcout << "starting simulation\n"; force_output();

  for (int t = 0; t <= maxTime; ++t) {
    if (is_in_time_points(t, time_points)) {
   //   Rcout << "update_unphased\n"; force_output();
      O.update_unphased(Pop, t, record_true_junctions, numRecombinations,
                        num_indiv_sampled);
    }

    std::vector< Fish_inf > newGeneration(popSize);

    update_pop(Pop, newGeneration, popSize, numRecombinations, num_threads);

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
List sim_phased_unphased_cpp(int pop_size,
                             double freq_ancestor_1,
                             int total_runtime,
                             double size_in_morgan,
                             NumericVector markers,
                             NumericVector time_points,
                             int seed,
                             bool verbose,
                             bool record_true_junctions,
                             int num_indiv_sampled,
                             int num_threads) {

  rnd_t rndgen(seed);

  std::vector< double > marker_dist(markers.begin(), markers.end());

  Output O = simulation_phased_nonphased(pop_size, freq_ancestor_1,
                                         total_runtime,
                                         size_in_morgan,
                                         marker_dist,
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

int get_seed() {
  const auto tt = static_cast<int64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  auto tid = std::this_thread::get_id();
  const uint64_t e3{ std::hash<std::remove_const_t<decltype(tid)>>()(tid) };
  return static_cast<int>(tt + e3);
}