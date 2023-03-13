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
};

using chrom = std::vector<junction>;


struct Fish_inf {
    chrom chromosome1;
    chrom chromosome2;

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



struct Fish_multi {
    std::vector< chrom > genome1;
    std::vector< chrom > genome2;

    Fish_multi(const int initLoc, int num_chrom) {
        junction left = junction(0.0, initLoc);
        junction right = junction(1, -1);
        chrom local_chrom = {left ,right};

        //genome1 = std::vector< chrom >(num_chrom, local_chrom);
        //genome2 = genome1;

        for (size_t i = 0; i < num_chrom; ++i) {
            genome1.push_back(local_chrom);
            genome2.push_back(local_chrom);
        }
    }

    Fish_multi(int num_chrom) {
        genome1 = std::vector< chrom >(num_chrom);
        genome2 = std::vector< chrom >(num_chrom);
    }

    Fish_multi(const Fish_multi& other) {
        genome1 = other.genome1;
        genome2 = other.genome2;
    }
};

Fish_fin mate_fin(const Fish_fin& A, const Fish_fin& B,
                  double numRecombinations, rnd_t& rndgen);

Fish_inf mate(const Fish_inf& A, const Fish_inf& B,
                  const std::vector<double>& numRecombinations, rnd_t& rndgen);

Fish_multi mate(const Fish_multi& A, const Fish_multi& N,
                const std::vector<double>& numRecombinations, rnd_t& rndgen);

long double getRecomPos();
int getRecomPos(int L, rnd_t& rndgen);

bool is_in_time_points(int t,
                       const Rcpp::NumericVector& time_points);



#endif /* Fish_hpp */
