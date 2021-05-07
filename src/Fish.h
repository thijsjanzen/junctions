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

struct Fish_explicit {
    std::vector< int > chromosome1;
    std::vector< int > chromosome2;

    Fish_explicit()
    {}

    Fish_explicit(const std::vector< int >& c1,
                  const std::vector< int >& c2) :
        chromosome1(c1), chromosome2(c2) {
    }

    Fish_explicit(int allele, size_t num_markers) {
        chromosome1 = std::vector<int>(num_markers, allele);
        chromosome2 = chromosome1;
    }

    std::vector< int > gamete(double morgan,
                              rnd_t& rndgen,
                              const emp_genome& emp_gen) const {

        std::vector<size_t> recom_pos = emp_gen.recompos(morgan,
                                                         rndgen);

        if (recom_pos.size() == 1) {
            if(rndgen.random_number(2)) {
                return chromosome1;
            }
            return chromosome2;
        }

        std::vector < std::vector<int>::const_iterator > iters = {chromosome1.begin(),
                                                                  chromosome2.begin()};
        std::vector< int > recombined_chromosome;
        int index = rndgen.random_number(2);
        size_t prev_start = 0;

        for(size_t i = 0; i < recom_pos.size(); ++i) {
            auto start = iters[index] + prev_start;
            auto end   = iters[index] + recom_pos[i];

            prev_start = recom_pos[i];
            recombined_chromosome.insert(recombined_chromosome.end(), start, end);
            index = 1 - index;
        }

        return recombined_chromosome;
    }
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
