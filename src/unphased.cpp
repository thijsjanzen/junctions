#include <stdio.h>
#include "Output.h"
#include "Fish.h"
#include "random_functions.h"

#include <Rcpp.h>
using namespace Rcpp;

Output simulation_phased_nonphased(int popSize,
                                   double initRatio,
                                   int maxTime,
                                   double numRecombinations,
                                   int numberOfMarkers,
                                   const NumericVector& time_points,
                                   bool verbose
)    {

  Output O;
  std::vector< Fish_inf > Pop;
  std::vector<double> markers = generate_random_markers(numberOfMarkers);

  O.markers = markers;

  Fish_inf parent1 = Fish_inf(0);
  Fish_inf parent2 = Fish_inf(1);

  for(int i = 0; i < popSize; ++i) {
    Fish_inf p1 = parent1;
    Fish_inf p2 = parent1;

    if(uniform() < initRatio) {
      p1 = parent2;
    }
    if(uniform() < initRatio) {
      p2 = parent2;
    }

    Pop.push_back( mate_inf(p1,p2, numRecombinations));
  }

  if(verbose) Rcout << "0--------25--------50--------75--------100\n";
  if(verbose) Rcout << "*";
  int updateFreq = maxTime / 20;
  if(updateFreq < 1) updateFreq = 1;

  for(int t = 0; t <= maxTime; ++t) {
    if(is_in_time_points(t, time_points)) {
      O.update_unphased(Pop, t);
    }

    std::vector< Fish_inf > newGeneration;

    for(int i = 0; i < popSize; ++i)  {
      int index1 = random_number(popSize);
      int index2 = random_number(popSize);
      while(index2 == index1) index2 = random_number(popSize);

      Fish_inf kid = mate_inf(Pop[index1], Pop[index2], numRecombinations);

      newGeneration.push_back(kid);
    }

    Pop = newGeneration;
    newGeneration.clear();
    if(verbose) {
      if(t % updateFreq == 0) {
        Rcout << "**";
      }
    }

    Rcpp::checkUserInterrupt();
  }
  if(verbose) Rcout << "\n";
  return O;
}

// [[Rcpp::export]]
List sim_phased_unphased_cpp(int pop_size,
                             double freq_ancestor_1,
                             int total_runtime,
                             double size_in_morgan,
                             int number_of_markers,
                             NumericVector time_points,
                             int seed,
                             bool verbose) {

  set_seed(seed);

  Output O = simulation_phased_nonphased(pop_size, freq_ancestor_1,
                                         total_runtime,
                                         size_in_morgan, number_of_markers,
                                         time_points,
                                         verbose);

  int num_rows = O.results.size();
  int num_cols = O.results[0].size();

  NumericMatrix output_matrix(num_rows, num_cols);
  for(int i = 0; i < num_rows; ++i) {
    for(int j = 0; j < num_cols; ++j) {
      output_matrix(i,j) = O.results[i][j];
    }
  }

  return List::create(Named("results") = output_matrix);
}