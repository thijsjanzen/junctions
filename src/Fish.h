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
//  Fish.hpp
//
//
//  Created by Thijs Janzen on 07/11/2017.
//
//

#pragma once

#include <stdio.h>
#include <vector>
#include <algorithm>
#include "random_functions.h"   // NOLINT [build/include_subdir]
#include <Rcpp.h>

struct junction {
    double pos;
    int right;

    junction();
    junction(double loc, int A);
    junction(const junction& other);
};


struct Fish_inf {
    std::vector< junction > chromosome1;
    std::vector< junction > chromosome2;

    Fish_inf();
    explicit Fish_inf(int initLoc);
    Fish_inf(const Fish_inf& other);
    Fish_inf(Fish_inf&& other);
    Fish_inf& operator=(Fish_inf&& other);
    Fish_inf& operator=(const Fish_inf& other);
};

struct Fish_fin  {
    std::vector<bool> chromosome1;
    std::vector<bool> chromosome2;

    Fish_fin();
    Fish_fin(const bool initLoc, const int genomeSize);
};

Fish_fin mate_fin(const Fish_fin& A, const Fish_fin& B,
                  double numRecombinations, rnd_t* rndgen);

Fish_inf mate_inf(const Fish_inf& A, const Fish_inf& B,
                  double numRecombinations, rnd_t* rndgen);

long double getRecomPos();
int getRecomPos(int L, rnd_t* rndgen);

bool is_in_time_points(int t,
                       const Rcpp::NumericVector& time_points);
