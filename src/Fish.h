//
//  Fish.hpp
//
//
//  Created by Thijs Janzen on 07/11/2017.
//
//

#ifndef Fish_hpp
#define Fish_hpp

#include <stdio.h>
#include <vector>
#include <algorithm>
#include "random_functions.h"
#include <Rcpp.h>

struct junction {
    double pos;
    int right;

    junction();
    junction(double loc, int A);
    junction(const junction& other);

    bool operator <(const junction& other) const;
};


struct Fish_inf {
    std::vector< junction > chromosome1;
    std::vector< junction > chromosome2;

    Fish_inf();
    Fish_inf(int initLoc);
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
          double numRecombinations, rnd_t& rndgen);

Fish_inf mate_inf(const Fish_inf& A, const Fish_inf& B,
                  double numRecombinations, rnd_t& rndgen);

long double getRecomPos();
int getRecomPos(int L, rnd_t& rndgen);

bool is_in_time_points(int t,
                       const Rcpp::NumericVector& time_points);

#endif /* Fish_hpp */
