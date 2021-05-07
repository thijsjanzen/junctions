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

void update_pop(const std::vector<Fish_explicit>& old_pop,
                std::vector<Fish_explicit>& pop,
                int popSize,
                int morgan,
                int num_threads,
                const emp_genome& emp_gen) {

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

        pop[i] = Fish_explicit(old_pop[index1].gamete(morgan, rndgen2, emp_gen),
                              old_pop[index2].gamete(morgan, rndgen2, emp_gen));
      }
    }
  );
}


Output simulation_phased_nonphased(int popSize,
                                   double initRatio,
                                   int maxTime,
                                   double numRecombinations,
                                   const std::vector< double >&  markers,
                                   const NumericVector& time_points,
                                   bool verbose,
                                   bool record_true_junctions,
                                   int num_indiv_sampled,
                                   int num_threads,
                                   rnd_t& rndgen)    {

  Output O;
  std::vector< Fish_explicit > Pop(popSize);

  O.markers = markers;
  emp_genome empgen(markers);

  Fish_explicit parent1 = Fish_explicit(0, markers.size());
  Fish_explicit parent2 = Fish_explicit(1, markers.size());

  for (int i = 0; i < popSize; ++i) {
    Fish_explicit p1 = parent2;
    Fish_explicit p2 = parent2;

    if (rndgen.uniform() < initRatio) {
      p1 = parent1;
    }
    if (rndgen.uniform() < initRatio) {
      p2 = parent1;
    }

    Pop[i] = Fish_explicit(p1.gamete(numRecombinations, rndgen, empgen),
                           p2.gamete(numRecombinations, rndgen, empgen)) ;   ///mate_inf(p1,p2, numRecombinations, rndgen);
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

    std::vector< Fish_explicit > newGeneration(popSize);

    update_pop(Pop, newGeneration, popSize, numRecombinations, num_threads,
               empgen);

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
List sim_phased_unphased_explicit_cpp(int pop_size,
                             double freq_ancestor_1,
                             int total_runtime,
                             double size_in_morgan,
                             NumericVector markers,
                             NumericVector time_points,
                             bool verbose,
                             bool record_true_junctions,
                             int num_indiv_sampled,
                             int num_threads) {

  rnd_t rndgen;

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

  return List::create(Named("results") = output_matrix);
}