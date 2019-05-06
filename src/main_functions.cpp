//
//  main_functions.cpp
//
//
//  Created by Thijs Janzen on 07/11/2017.
//
//
#include <stdio.h>
#include "Output.h"
#include "Fish.h"
#include "random_functions.h"

#include <Rcpp.h>
using namespace Rcpp;

Output doSimulation_inf(int popSize,
                    double initRatio,
                    int maxTime,
                    double numRecombinations,
                    int numberOfMarkers)    {
    Output O;
    std::vector<Fish_inf> Pop;
    std::vector<double> markers;
    if(numberOfMarkers > 0) {
        for(int i = 0; i < numberOfMarkers; ) {
            double pos = uniform();
            if(pos > 0 && pos < 1.0) {
                ++i;
                markers.push_back(pos);
            }
        }
        std::sort(markers.begin(), markers.end());
    }

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

        Pop.push_back(mate_inf(p1,p2, numRecombinations));
    }

    int updateFreq = maxTime / 20;
    if(updateFreq < 1) updateFreq = 1;

    for(int t = 0; t < maxTime; ++t) {
        O.update_inf(Pop);
        if(numberOfMarkers > 0) O.detectNumJunctions(Pop, markers);

        std::vector<Fish_inf> newGeneration;

        for(int i = 0; i < popSize; ++i)  {
            int index1 = random_number(popSize);
            int index2 = random_number(popSize);

            Fish_inf kid = mate_inf(Pop[index1], Pop[index2], numRecombinations);

            newGeneration.push_back(kid);
        }

        Pop = newGeneration;
        newGeneration.clear();
    }
    return O;
}

Output doSimulation_fin(int popSize,
                        int genomeSize,
                        double initRatio,
                        int maxTime,
                        double numRecombinations)    {
    Output O;
    std::vector<Fish_fin> Pop;

    Fish_fin parent1 = Fish_fin(0, genomeSize);
    Fish_fin parent2 = Fish_fin(1, genomeSize);

    for(int i = 0; i < popSize; ++i) {
        Fish_fin p1 = parent1;
        Fish_fin p2 = parent1;

        if(uniform() < initRatio) {
            p1 = parent2;
        }
        if(uniform() < initRatio) {
            p2 = parent2;
        }

        Pop.push_back(mate_fin(p1,p2, numRecombinations));
    }

    int updateFreq = maxTime / 20;
    if(updateFreq < 1) updateFreq = 1;

    for(int t = 0; t < maxTime; ++t) {
        O.update_fin(Pop);
        std::vector<Fish_fin> newGeneration;

        for(int i = 0; i < popSize; ++i)
        {
            int index1 = random_number(popSize);
            int index2 = random_number(popSize);

            Fish_fin kid;
            kid = mate_fin(Pop[index1],Pop[index2], numRecombinations);

            newGeneration.push_back(kid);
        }

        Pop = newGeneration;
        newGeneration.clear();
    }

    return O;
}

// [[Rcpp::export]]
List sim_fin_chrom(int pop_size,
                   double freq_ancestor_1,
                   int run_time,
                   double size_in_Morgan,
                   int seed,
                   int R) {

  set_seed(seed);

  Output O = doSimulation_fin(pop_size,
                              R + 1,
                              freq_ancestor_1,
                              run_time,
                              size_in_Morgan);

  return List::create(Named("avgJunctions") = O.avgJunctions);
}

// [[Rcpp::export]]
List sim_inf_chrom(int pop_size,
                   double freq_ancestor_1,
                   int run_time,
                   double size_in_Morgan,
                   int markers,
                   int seed) {
    set_seed(seed);

    //Rcout << "sim_inf_chrom, let's go!\n";

    Output O = doSimulation_inf(pop_size,
                                freq_ancestor_1,
                                run_time,
                                size_in_Morgan,
                                markers);
    return List::create(Named("avgJunctions") = O.avgJunctions,
                        Named("detectedJunctions") = O.avg_detected_Junctions,
                        Named("markers") = O.markers);
}
