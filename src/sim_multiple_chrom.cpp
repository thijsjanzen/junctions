#include <stdio.h>
#include <vector>
#include "Output.h"
#include "Fish.h"
#include "random_functions.h"
#include "util.h"

#include <thread>
#include <chrono>
#include <functional>

#include <RcppParallel.h>
#include <Rcpp.h>
using namespace Rcpp;


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

  O.markers_multi = markers;

  Fish_multi parent1 = Fish_multi(0, numRecombinations.size());
  Fish_multi parent2 = Fish_multi(1, numRecombinations.size());

  std::vector< Fish_multi > Pop(popSize, parent1);

  for (int i = 0; i < popSize; ++i) {
    Fish_multi p1 = parent2;
    Fish_multi p2 = parent2;

    if (rndgen.uniform() < initRatio) {
      p1 = parent1;
    }
    if (rndgen.uniform() < initRatio) {
      p2 = parent1;
    }

    Pop[i] = mate(p1, p2, numRecombinations, rndgen);
  }

  if (verbose) Rcout << "0--------25--------50--------75--------100\n";
  if (verbose) Rcout << "*";
  int updateFreq = maxTime / 20;
  if (updateFreq < 1) updateFreq = 1;

//  Rcout << "starting simulation\n"; force_output();
  std::vector< Fish_multi > newGeneration = Pop;

  for (size_t t = 0; t <= maxTime; ++t) {
    if (is_in_time_points(t, time_points)) {
      O.update_unphased(Pop, t, record_true_junctions, numRecombinations,
                        num_indiv_sampled);
    }

    update_pop<Fish_multi>(Pop, newGeneration, popSize, numRecombinations, num_threads);

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
                                   const Rcpp::NumericMatrix& markers,
                                   const Rcpp::NumericVector& time_points,
                                   bool verbose,
                                   bool record_true_junctions,
                                   int num_indiv_sampled,
                                   int num_threads) {

  rnd_t rndgen;

 // std::vector< double > marker_dist(markers.begin(), markers.end());

  std::vector< std::vector< double >> markers_cpp(markers.nrow(),
                                                  std::vector<double>(markers.ncol(), 0.0));
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

  int num_rows = O.results_multi.size();
  int num_cols = O.results_multi[0].size();

  NumericMatrix output_matrix(num_rows, num_cols);
  for (int i = 0; i < num_rows; ++i) {
    for (int j = 0; j < num_cols; ++j) {
      output_matrix(i,j) = O.results_multi[i][j];
    }
  }

  if (record_true_junctions) {
    int num_rows_t = O.true_results_multi.size();
    int num_cols_t = O.true_results_multi[0].size();
    NumericMatrix output_matrix_true(num_rows_t, num_cols_t);
    for(int i = 0; i < num_rows_t; ++i) {
      for(int j = 0; j < num_cols_t; ++j) {
        output_matrix_true(i, j) = O.true_results_multi[i][j];
      }
    }

    return List::create(Named("results") = output_matrix,
                        Named("true_results") = output_matrix_true);
  }
  return List::create(Named("results") = output_matrix);
}