//
//  Output.hpp
//
//
//  Created by Thijs Janzen on 07/11/2017.
//
//

#ifndef Output_hpp
#define Output_hpp

#include <stdio.h>
#include "Fish.h"
#include <vector>

#include <chrono>
#include <thread>


struct Output {
    std::vector<double> avgJunctions;
    std::vector<double> avg_detected_Junctions;
    std::vector<double> markers;
    // the average heterozygosity at t
    std::vector<double> avg_hetero;

    // the full distribution of junctions, at t
    std::vector< std::vector< int > > junction_dist;

    std::vector< std::vector< double > > results;
    std::vector< std::vector< double > > true_results;

    void update_inf(const std::vector< Fish_inf >& Pop);
    void update_fin(const std::vector< Fish_fin >& Pop);

    void update_unphased(const std::vector< Fish_inf >& Pop,
                         size_t t,
                         bool record_true_junctions,
                         double morgan,
                         size_t num_indiv);

    void update_unphased(const std::vector< Fish_explicit >& Pop,
                                 size_t t,
                                 bool record_true_junctions,
                                 double morgan,
                                 size_t num_indiv);


    void detectNumJunctions(const std::vector<Fish_inf> &Pop,
                            const std::vector<double> &markers);

    void detect_junctions_backcross(const std::vector< Fish_inf > &Pop,
                                      const std::vector<double> &markers);
};

#endif /* Output_hpp */
