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
//  Output.hpp
//
//
//  Created by Thijs Janzen on 07/11/2017.
//
//

#pragma once

#include <stdio.h>
#include "Fish.h"
#include <vector>

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
