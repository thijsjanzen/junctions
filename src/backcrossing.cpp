#include <stdio.h>
#include "Output.h"
#include "Fish.h"
#include "random_functions.h"

#include <Rcpp.h>
using namespace Rcpp;

Output doSimulation_backcrossing(int population_size,
                                 int maximum_time,
                                 double number_of_recombinations,
                                 int number_of_markers,
                                 const NumericVector& time_points)    {

  // declaration of data holders
  Output O;
  std::vector< Fish_inf > Pop;
  std::vector<double> markers;

  // generate random markers if necessary
  if(number_of_markers > 0) {
    markers = generate_random_markers(number_of_markers);
  }
  O.markers = markers;

  Fish_inf parent1 = Fish_inf(0);
  Fish_inf parent2 = Fish_inf(1);

  Fish_inf back_cross_parent = parent1;

  // initialize the simulation with individuals that are all a cross of parent1 x parent2
  // in principle these are all identical, so this is a bit inefficient coding,
  // but shouldn't cost much time anyway.
  for(int i = 0; i < population_size; ++i) {
    Pop.push_back(mate_inf(parent1, parent2,
                           number_of_recombinations));
  }

  // because we initialize with F1, we start at t = 1
  for(int t = 0; t < maximum_time; ++t) {
    // record all information
    if(is_in_time_points(t, time_points)) {
      O.update_inf(Pop);
      O.detect_junctions_backcross(Pop, markers);
    }

    // generate the next generation. We assume non-overlapping generations
    std::vector< Fish_inf > next_generation;

    for(int i = 0; i < population_size; ++i)  {
      int index1 = random_number(population_size);

      Fish_inf kid = mate_inf(Pop[index1], back_cross_parent,
                              number_of_recombinations);

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
List simulate_backcrossing_cpp(int population_size,
                               int maximum_time,
                               double size_in_morgan,
                               int markers,
                               int seed,
                               NumericVector time_points) {

  set_seed(seed);

  Output O = doSimulation_backcrossing(population_size,
                                       maximum_time,
                                       size_in_morgan,
                                       markers,
                                       time_points);

  return List::create(Named("average_junctions") = O.avgJunctions,
                      Named("detected_junctions") = O.avg_detected_Junctions,
                      Named("markers") = O.markers,
                      Named("junction_distribution") = O.junction_dist,
                      Named("average_heterozygosity") = O.avg_hetero);
}

