// Copyright 2018 - 2024 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//

#include <stdio.h>
#include <vector>

#include "tbb_stub.h"          // NOLINT [build/include_subdir]

#include "Output.h"            // NOLINT [build/include_subdir]
#include "Fish.h"              // NOLINT [build/include_subdir]
#include "random_functions.h"  // NOLINT [build/include_subdir]

#include <RcppParallel.h>
#include <Rcpp.h>

// Parents and per-offspring seeds are drawn on the main thread from the
// seeded master generator; each offspring is then produced with its own
// generator. This makes results reproducible from the seed, independently of
// the number of threads and of how TBB splits the work.
void update_pop(const std::vector<Fish_inf>& old_pop,
                std::vector<Fish_inf>* pop,
                int popSize,
                double numRecombinations,
                rnd_t* rndgen,
                tbb::task_arena* arena) {
  std::vector<int> parent1(popSize);
  std::vector<int> parent2(popSize);
  std::vector<unsigned int> seeds(popSize);
  for (int i = 0; i < popSize; ++i) {
    parent1[i] = rndgen->random_number(popSize);
    parent2[i] = rndgen->random_number(popSize);
    while (parent2[i] == parent1[i]) {
      parent2[i] = rndgen->random_number(popSize);
    }
    seeds[i] = static_cast<unsigned int>(rndgen->rndgen_());
  }
  arena->execute([&] {
    tbb::parallel_for<int>(0, popSize, 1, [&](int i) {
      rnd_t local_rndgen(seeds[i]);
      (*pop)[i] = mate_inf(old_pop[parent1[i]], old_pop[parent2[i]],
                           numRecombinations, &local_rndgen);
    });
  });
}

Output simulation_phased_nonphased(int popSize,
                                   double initRatio,
                                   int maxTime,
                                   double numRecombinations,
                                   std::vector< double >  markers,
                                   const Rcpp::NumericVector& time_points,
                                   bool verbose,
                                   bool record_true_junctions,
                                   int num_indiv_sampled,
                                   rnd_t* rndgen,
                                   int num_threads)    {
  Output O;
  std::vector< Fish_inf > Pop(popSize);
  tbb::task_arena arena(num_threads);

  O.markers = markers;

  Fish_inf parent1 = Fish_inf(0);
  Fish_inf parent2 = Fish_inf(1);

  for (int i = 0; i < popSize; ++i) {
    Fish_inf p1 = parent2;
    Fish_inf p2 = parent2;

    if (rndgen->uniform() < initRatio) {
      p1 = parent1;
    }
    if (rndgen->uniform() < initRatio) {
      p2 = parent1;
    }

    Pop[i] = mate_inf(p1, p2, numRecombinations, rndgen);
  }

  if (verbose) Rcpp::Rcout << "0--------25--------50--------75--------100\n";
  if (verbose) Rcpp::Rcout << "*";
  int updateFreq = maxTime / 20;
  if (updateFreq < 1) updateFreq = 1;

  for (int t = 0; t <= maxTime; ++t) {
    if (is_in_time_points(t, time_points)) {
      O.update_unphased(Pop, t, record_true_junctions, numRecombinations,
                        num_indiv_sampled);
    }

    std::vector< Fish_inf > newGeneration(popSize);

    update_pop(Pop, &newGeneration, popSize, numRecombinations, rndgen,
               &arena);

    Pop.swap(newGeneration);

    if (verbose) {
      if (t % updateFreq == 0) {
        Rcpp::Rcout << "**";
      }
    }
    Rcpp::checkUserInterrupt();
  }
  if (verbose) Rcpp::Rcout << "\n";
  return O;
}

// [[Rcpp::export]]
Rcpp::List sim_phased_unphased_cpp(int pop_size,
                                   double freq_ancestor_1,
                                   int total_runtime,
                                   double size_in_morgan,
                                   Rcpp::NumericVector markers,
                                   Rcpp::NumericVector time_points,
                                   bool verbose,
                                   bool record_true_junctions,
                                   int num_indiv_sampled,
                                   int seed) {
  int num_threads = get_rcpp_num_threads();

  rnd_t rndgen(static_cast<unsigned int>(seed));
  std::vector< double > marker_dist(markers.begin(), markers.end());
  Output O = simulation_phased_nonphased(pop_size,
                                         freq_ancestor_1,
                                         total_runtime,
                                         size_in_morgan,
                                         marker_dist,
                                         time_points,
                                         verbose,
                                         record_true_junctions,
                                         num_indiv_sampled,
                                         &rndgen,
                                         num_threads);
  if (O.results.empty()) {
    Rcpp::stop("no data recorded: none of time_points falls within "
               "0 .. total_runtime");
  }
  int num_rows = O.results.size();
  int num_cols = O.results[0].size();

  Rcpp::NumericMatrix output_matrix(num_rows, num_cols);
  for (int i = 0; i < num_rows; ++i) {
    for (int j = 0; j < num_cols; ++j) {
      output_matrix(i, j) = O.results[i][j];
    }
  }

  if (record_true_junctions && !O.true_results.empty()) {
    int num_rows_t = O.true_results.size();
    int num_cols_t = O.true_results[0].size();
    Rcpp::NumericMatrix output_matrix_true(num_rows_t, num_cols_t);
    for (int i = 0; i < num_rows_t; ++i) {
      for (int j = 0; j < num_cols_t; ++j) {
        output_matrix_true(i, j) = O.true_results[i][j];
      }
    }

    return Rcpp::List::create(Rcpp::Named("results") = output_matrix,
                              Rcpp::Named("true_results") = output_matrix_true);
  }
  return Rcpp::List::create(Rcpp::Named("results") = output_matrix);
}
