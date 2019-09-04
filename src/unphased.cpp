#include <stdio.h>
#include "Output.h"
#include "Fish.h"
#include "random_functions.h"

#include <omp.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]


Output simulation_phased_nonphased_threaded(int popSize,
                                            double initRatio,
                                            int maxTime,
                                            double numRecombinations,
                                            int numberOfMarkers,
                                            const NumericVector& time_points,
                                            bool verbose,
                                            int num_threads)    {

  Output O;
  std::vector< Fish_inf > Pop;
  std::vector<double> markers = generate_random_markers(numberOfMarkers);

  O.markers = markers;

  Fish_inf parent1 = Fish_inf(0);
  Fish_inf parent2 = Fish_inf(1);

  for(int i = 0; i < popSize; ++i) {
    Fish_inf p1 = parent2;
    Fish_inf p2 = parent2;

    if(uniform() < initRatio) {
      p1 = parent1;
    }
    if(uniform() < initRatio) {
      p2 = parent1;
    }

    Pop.push_back( mate_inf(p1,p2, numRecombinations));
  }

  omp_set_num_threads(num_threads);

//  int established_num_threads = omp_get_num_threads();
 // Rcout << "running with " << established_num_threads << " threads\n";


  if(verbose) Rcout << "0--------25--------50--------75--------100\n";
  if(verbose) Rcout << "*";
  int updateFreq = maxTime / 20;
  if(updateFreq < 1) updateFreq = 1;



  for(int t = 0; t <= maxTime; ++t) {
    if(is_in_time_points(t, time_points)) {
      O.update_unphased(Pop, t);
    }

    std::vector< Fish_inf > newGeneration(popSize);

#pragma omp parallel for shared(newGeneration, popSize, Pop)
    for(size_t i = 0; i < popSize; ++i)  {

// #pragma omp critical
// {
//       Rcout << t << " " << omp_get_thread_num() << " " << i << " " << omp_get_num_threads() << "\n";
// }
      int index1 = random_number(popSize);
      int index2 = random_number(popSize);
      while(index2 == index1) index2 = random_number(popSize);

      Fish_inf kid = mate_inf(Pop[index1], Pop[index2], numRecombinations);

      newGeneration[i] = kid;
    }

    Pop.swap(newGeneration);

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


//' simulate junctions
//' @param pop_size Population Size
//' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
//' @param total_runtime Maximum time after which the simulation is to be stopped
//' @param size_in_morgan Mean number of crossovers per meiosis (e.g. size in Morgan of the chromosome)
//' @param number_of_markers The number of genetic markers superimposed on the chromosome.
//' @param time_points vector with time points at which local ancestry has to be recorded to be returned at the end of the simulation. If left at -1, ancestry is recorded at every generation (computationally heavy).
//' @param seed Seed of the pseudo-random number generator
//' @param verbose displays a progress bar
//' @param num_threads if larger than one, multithreading is used.
//' @export
// [[Rcpp::export]]
List sim_phased_unphased_threaded_cpp(int pop_size,
                                      double freq_ancestor_1,
                                      int total_runtime,
                                      double size_in_morgan,
                                      int number_of_markers,
                                      NumericVector time_points,
                                      int seed,
                                      bool verbose,
                                      int num_threads) {

  set_seed(seed);

  Output O = simulation_phased_nonphased_threaded(pop_size, freq_ancestor_1,
                                                  total_runtime,
                                                  size_in_morgan, number_of_markers,
                                                  time_points,
                                                  verbose,
                                                  num_threads);

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



Output simulation_phased_nonphased(int popSize,
                                   double initRatio,
                                   int maxTime,
                                   double numRecombinations,
                                   int numberOfMarkers,
                                   const NumericVector& time_points,
                                   bool verbose
)    {

  Output O;
  std::vector< Fish_inf > Pop(popSize);
  std::vector<double> markers = generate_random_markers(numberOfMarkers);

  O.markers = markers;

  Fish_inf parent1 = Fish_inf(0);
  Fish_inf parent2 = Fish_inf(1);

  for(int i = 0; i < popSize; ++i) {
    Fish_inf p1 = parent2;
    Fish_inf p2 = parent2;

    if(uniform() < initRatio) {
      p1 = parent1;
    }
    if(uniform() < initRatio) {
      p2 = parent1;
    }

    Pop[i] = mate_inf(p1,p2, numRecombinations);
  }

  if(verbose) Rcout << "0--------25--------50--------75--------100\n";
  if(verbose) Rcout << "*";
  int updateFreq = maxTime / 20;
  if(updateFreq < 1) updateFreq = 1;

  for(int t = 0; t <= maxTime; ++t) {
    if(is_in_time_points(t, time_points)) {
      O.update_unphased(Pop, t);
    }

    std::vector< Fish_inf > newGeneration(popSize);

    for(int i = 0; i < popSize; ++i)  {
      int index1 = random_number(popSize);
      int index2 = random_number(popSize);
      while(index2 == index1) index2 = random_number(popSize);

      Fish_inf kid = mate_inf(Pop[index1], Pop[index2], numRecombinations);

      newGeneration[i] = kid;
    }

    Pop.swap(newGeneration);

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



//' simulate junctions
//' @param pop_size Population Size
//' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
//' @param total_runtime Maximum time after which the simulation is to be stopped
//' @param size_in_morgan Mean number of crossovers per meiosis (e.g. size in Morgan of the chromosome)
//' @param number_of_markers The number of genetic markers superimposed on the chromosome.
//' @param time_points vector with time points at which local ancestry has to be recorded to be returned at the end of the simulation. If left at -1, ancestry is recorded at every generation (computationally heavy).
//' @param seed Seed of the pseudo-random number generator
//' @param verbose displays a progress bar
//' @export
//' @export
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
