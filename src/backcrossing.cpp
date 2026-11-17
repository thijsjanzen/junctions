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

#include "Output.h"             // NOLINT [build/include_subdir]
#include "Fish.h"               // NOLINT [build/include_subdir]
#include "random_functions.h"   // NOLINT [build/include_subdir]

#include <Rcpp.h>

Output doSimulation_backcrossing(int population_size,
                                 double freq_ancestor_1,
                                 int total_runtime,
                                 double size_in_morgan,
                                 int number_of_markers,
                                 const Rcpp::NumericVector& time_points,
                                 rnd_t* rndgen)    {
  // declaration of data holders
  Output O;
  std::vector< Fish_inf > Pop;
  std::vector<double> markers;

  // generate random markers if necessary
  if (number_of_markers > 0) {
    markers = rndgen->generate_random_markers(number_of_markers);
  }
  O.markers = markers;

  Fish_inf parent1 = Fish_inf(0);
  Fish_inf parent2 = Fish_inf(1);

  Fish_inf back_cross_parent = parent1;

  // initialize the simulation with individuals that are all
  // a cross of parent1 x parent2
  // in principle these are all identical, so this is a bit inefficient coding,
  // but shouldn't cost much time anyway.
  for (int i = 0; i < population_size; ++i) {
    Fish_inf focal_parent1 = parent1;
    Fish_inf focal_parent2 = parent2;

    Pop.push_back(mate_inf(focal_parent1, focal_parent2,
                           size_in_morgan, rndgen));
  }

  // because we initialize with F1, we start at t = 1
  for (int t = 0; t < total_runtime; ++t) {
    // record all information
    if (is_in_time_points(t, time_points)) {
      O.update_inf(Pop);
      O.detect_junctions_backcross(Pop, markers);
    }

    // generate the next generation. We assume non-overlapping generations
    std::vector< Fish_inf > next_generation;

    for (int i = 0; i < population_size; ++i)  {
      int index1 = rndgen->random_number(population_size);

      Fish_inf kid = mate_inf(Pop[index1], back_cross_parent,
                              size_in_morgan, rndgen);

      next_generation.push_back(kid);
    }

    Pop = next_generation;
    next_generation.clear();  // just to be sure

    // checkUserInterrupt allows for the user to hit the 'stop' button in R
    Rcpp::checkUserInterrupt();
  }

  // return the recorded data
  return O;
}

// [[Rcpp::export]]
Rcpp::List simulate_backcrossing_cpp(int pop_size,
                               double freq_ancestor_1,
                               int total_runtime,
                               double size_in_morgan,
                               int number_of_markers,
                               Rcpp::NumericVector time_points,
                               int seed) {
  rnd_t rndgen;
  rndgen.set_seed(seed);

  Output O = doSimulation_backcrossing(pop_size,
                                       freq_ancestor_1,
                                       total_runtime,
                                       size_in_morgan,
                                       number_of_markers,
                                       time_points,
                                       &rndgen);

  return Rcpp::List::create(Rcpp::Named("average_junctions") =
                              O.avgJunctions,
                            Rcpp::Named("detected_junctions") =
                              O.avg_detected_Junctions,
                            Rcpp::Named("markers") =
                              O.markers,
                            Rcpp::Named("junction_distribution") =
                              O.junction_dist,
                            Rcpp::Named("average_heterozygosity") =
                              O.avg_hetero);
}

