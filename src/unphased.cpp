#include <stdio.h>
#include <vector>
#include "Output.h"
#include "Fish.h"
#include "random_functions.h"

#include <tbb/tbb.h>

#include <Rcpp.h>
using namespace Rcpp;

Output simulation_phased_nonphased(int popSize,
                                   double initRatio,
                                   int maxTime,
                                   double numRecombinations,
                                   std::vector< double >  markers,
                                   const NumericVector& time_points,
                                   bool verbose,
                                   bool record_true_junctions,
                                   int num_indiv_sampled,
                                   int num_threads)    {

  Output O;
  std::vector< Fish_inf > Pop(popSize);

  O.markers = markers;

  Fish_inf parent1 = Fish_inf(0);
  Fish_inf parent2 = Fish_inf(1);

  tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);

  for (int i = 0; i < popSize; ++i) {
    Fish_inf p1 = parent2;
    Fish_inf p2 = parent2;

    if (uniform() < initRatio) {
      p1 = parent1;
    }
    if (uniform() < initRatio) {
      p2 = parent1;
    }

    Pop[i] = mate_inf(p1,p2, numRecombinations);
  }

  if (verbose) Rcout << "0--------25--------50--------75--------100\n";
  if (verbose) Rcout << "*";
  int updateFreq = maxTime / 20;
  if (updateFreq < 1) updateFreq = 1;

  for (int t = 0; t <= maxTime; ++t) {
    if (is_in_time_points(t, time_points)) {
      O.update_unphased(Pop, t, record_true_junctions, numRecombinations,
                        num_indiv_sampled);
    }

    std::vector< Fish_inf > newGeneration(popSize);


    tbb::parallel_for(
      tbb::blocked_range<unsigned>(0, popSize),
      [&](const tbb::blocked_range<unsigned>& r) {
        for (unsigned i = r.begin(); i < r.end(); ++i) {

          int index1 = random_number(popSize);
          int index2 = random_number(popSize);
          while(index2 == index1) index2 = random_number_popsize();

          Fish_inf kid = mate_inf(Pop[index1], Pop[index2], numRecombinations);

          newGeneration[i] = kid;
        }
      }
    );

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

  set_seed(seed);
  set_poisson(size_in_morgan);
  set_random_number_popsize(pop_size);

  std::vector< double > marker_dist(markers.begin(), markers.end());

  Output O = simulation_phased_nonphased(pop_size, freq_ancestor_1,
                                         total_runtime,
                                         size_in_morgan,
                                         marker_dist,
                                         time_points,
                                         verbose,
                                         record_true_junctions,
                                         num_indiv_sampled,
                                         num_threads);

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
