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
//
//  Fish.cpp
//
//
//  Created by Thijs Janzen on 07/11/2017.
//
//

#include <utility>
#include <vector>
#include <algorithm>

#include "Fish.h"                           // NOLINT [build/include_subdir]
#include "random_functions.h"               // NOLINT [build/include_subdir]

#include <Rcpp.h>

int getRecomPos(int L,
                rnd_t* rnd) {
    int pos = -100;
    int index = rnd->random_number(L);
    // exclude the ends of the chromosome
    while (index == 0 || index == L) {
        index = rnd->random_number(L);
    }
    pos = index;

    return pos;
}

std::vector<junction> recombine_new(
    const std::vector<junction>& chromosome1,
    const std::vector<junction>& chromosome2,
    const std::vector<double>& recom_positions) {

    // we need something that is cheaply swappable:
    auto* g1 = &chromosome1;
    auto* g2 = &chromosome2;
    // offspring genome: recycle what's already there...
    std::vector<junction> go;


    // predicate for lower_bound
    auto less = [](const auto& j, double p) { return j.pos < p; };

    // helper lambda to get the value just *before* it.
    // we store the value to the right of a recombination-point
    // but we need the value to the left:
    auto value_at = [](auto begin, auto it) {
        return (begin != it) ? (it - 1)->right : -1;
    };

    double left_pos = 0.0;
    auto go_val = -1;
    for (auto right_pos : recom_positions) {
        auto it = std::lower_bound(g1->cbegin(), g1->cend(), left_pos, less);
        auto last = std::lower_bound(it, g1->cend(), right_pos, less);
        // [g1.first, it) : part of the genome *before* left_pos.
        // [it, last)     : part of the genome *after or equal to*
        //                  left_pos but *before* right_pos.
        auto g1_val = value_at(g1->cbegin(), it);
        if (g1_val != go_val) {
            if (it == last || it->pos != left_pos) {
                go.emplace_back(left_pos, g1_val);   // insert change to match
            } else {
                ++it;    // corner case: skip spurious double-change
            }
        }
        go.insert(go.end(), it, last);      // append [it, last)
        go_val = value_at(go.begin(), go.end());
        std::swap(g1, g2);
        left_pos = right_pos;
    }
    go.emplace_back(1.0, -1);
    return go;
}

std::vector<double> generate_recomPos(size_t number_of_recombinations,
                                      rnd_t* rndgen) {
    std::vector<double> recomPos(number_of_recombinations, 0);
    for (size_t i = 0; i < number_of_recombinations; ++i) {
        recomPos[i] = rndgen->uniform();
    }
    std::sort(recomPos.begin(), recomPos.end());
    if (recomPos.size() != number_of_recombinations) {
        Rcpp::stop("mismatch\n");
    }
    recomPos.push_back(1.0);

    return recomPos;
}

void Recombine_inf(std::vector<junction>* offspring,
                   const std::vector<junction>& chromosome1,
                   const std::vector<junction>& chromosome2,
                   double MORGAN,
                   rnd_t* rndgen)  {
    int numRecombinations = rndgen->poisson(MORGAN);

    if (numRecombinations == 0) {
        offspring->insert(offspring->end(),
                         chromosome1.begin(),
                         chromosome1.end());
        return;
    }
    std::vector<double> recomPos = generate_recomPos(numRecombinations,
                                                     rndgen);
    *offspring = recombine_new(chromosome1,
                              chromosome2,
                              recomPos);
    return;
}

Fish_inf mate_inf(const Fish_inf& A,
                  const Fish_inf& B,
                  double numRecombinations,
                  rnd_t* rndgen) {
    Fish_inf offspring;
    offspring.chromosome1.clear();
    offspring.chromosome2.clear();     // just to be sure.

    //  first the father chromosome
    int event = rndgen->random_number(2);
    switch (event) {
        case 0:  {
            Recombine_inf(&offspring.chromosome1,
                          A.chromosome1,
                          A.chromosome2,
                          numRecombinations,
                          rndgen);
            break;
        }
        case 1: {
            Recombine_inf(&offspring.chromosome1,
                          A.chromosome2,
                          A.chromosome1,
                          numRecombinations,
                          rndgen);
            break;
        }
    }

    //  then the mother chromosome
    event = rndgen->random_number(2);
    switch (event) {
        case 0:  {
            Recombine_inf(&offspring.chromosome2,
                          B.chromosome1,
                          B.chromosome2,
                          numRecombinations,
                          rndgen);
            break;
        }
        case 1: {
            Recombine_inf(&offspring.chromosome2,
                          B.chromosome2,
                          B.chromosome1,
                          numRecombinations,
                          rndgen);
            break;
        }
    }

    return offspring;
}

void Recombine_fin(std::vector<bool>* offspring,
                   std::vector<bool> chromosome1,
                   std::vector<bool> chromosome2,
                   double numberRecombinations,
                   rnd_t* rndgen)  {
    numberRecombinations = rndgen->poisson(numberRecombinations);

    // if there are not recombinations, preliminary exit
    if (numberRecombinations == 0) {
        offspring->insert(offspring->end(),
                          chromosome1.begin(),
                          chromosome1.end() );
        return;
    }

    std::vector<int> recomPos;
    // store L, so we avoid repeated calls of the function .size()
    int L = static_cast<int>(chromosome1.size());

    while (recomPos.size() < numberRecombinations) {
        int pos = getRecomPos(L, rndgen);
        recomPos.push_back(pos);
        // sort them, in case they are not sorted yet
        // we need this to remove duplicates, and later
        // to apply crossover
        std::sort(recomPos.begin(), recomPos.end());
        // remove duplicate recombination sites
        auto last = std::unique(recomPos.begin(), recomPos.end());
        recomPos.erase(last, recomPos.end());
    }

    // used to track which chromosome was used
    // during the last recombination event
    int order = 0;
    int start = 0;

    for (size_t i = 0; i < recomPos.size(); ++i) {
        int end = recomPos[i];
        if (order == 0) {  // add the first chromosome
            offspring->insert(offspring->end(),
                              chromosome1.begin() + start,
                              chromosome1.begin() + end);
            order = 1;
        } else {   // add the second chromosome
            offspring->insert(offspring->end(),
                              chromosome2.begin() + start,
                              chromosome2.begin() + end);
            order = 0;
        }
        start = end;
    }

    // add chromosomal content after
    // the last recombination site:
    if (order == 0) {
        offspring->insert(offspring->end(),
                          chromosome1.begin() + start,
                          chromosome1.end());
    } else {
        offspring->insert(offspring->end(),
                          chromosome2.begin() + start,
                          chromosome2.end());
    }

    return;
}

Fish_fin mate_fin(const Fish_fin& A,
                  const Fish_fin& B,
                  double numberRecombinations,
                  rnd_t* rndgen) {
    Fish_fin offspring;
    offspring.chromosome1.clear();
    offspring.chromosome2.clear();  // just to be sure.

    // random order or in other words,
    // we randomly select 1 of 2 produced chromosomes during recombination
    if (rndgen->uniform() < 0.5) {
        Recombine_fin(&offspring.chromosome1,
                      A.chromosome1, A.chromosome2,
                      numberRecombinations,
                      rndgen);
    } else {
        Recombine_fin(&offspring.chromosome1,
                      A.chromosome2, A.chromosome1,
                      numberRecombinations,
                      rndgen);
    }

    if (rndgen->uniform() < 0.5) {
        Recombine_fin(&offspring.chromosome2,
                      B.chromosome1, B.chromosome2,
                      numberRecombinations,
                      rndgen);
    } else {
        Recombine_fin(&offspring.chromosome2,
                      B.chromosome2, B.chromosome1,
                      numberRecombinations,
                      rndgen);
    }
    return offspring;
}
/////////////////////////////////////////////
//////////////  Member functions
/////////////////////////////////////////////

junction::junction() {
}

junction::junction(double loc, int A) : pos(loc), right(A) {
}

junction::junction(const junction& other) {
    pos = other.pos;
    right = other.right;
}

Fish_inf::Fish_inf() {
}

Fish_inf::Fish_inf(int initLoc)    {
    junction left  = junction(0.0, initLoc);
    junction right = junction(1, -1);
    chromosome1.push_back(left);
    chromosome1.push_back(right);
    chromosome2.push_back(left);
    chromosome2.push_back(right);
}

Fish_fin::Fish_fin() {
}

// constructor that sets all genome elements to "initLoc"
Fish_fin::Fish_fin(const bool initLoc, const int genomeSize) {
    chromosome1.clear();
    chromosome2.clear();
    chromosome1.resize(genomeSize, initLoc);
    chromosome2.resize(genomeSize, initLoc);
}

Fish_inf::Fish_inf(Fish_inf&& other) {
    chromosome1 = other.chromosome1;
    chromosome2 = other.chromosome2;
}

Fish_inf& Fish_inf::operator=(Fish_inf&& other) {
    if (this != &other) {
        chromosome1 = other.chromosome1;
        chromosome2 = other.chromosome2;
    }
    return *this;
}

Fish_inf::Fish_inf(const Fish_inf& other) {
    chromosome1 = other.chromosome1;
    chromosome2 = other.chromosome2;
}

Fish_inf& Fish_inf::operator=(const Fish_inf& other) {
    if (this != &other) {
        chromosome1 = other.chromosome1;
        chromosome2 = other.chromosome2;
    }
    return *this;
}


bool is_in_time_points(int t,
                       const Rcpp::NumericVector& time_points) {
    for (auto i : time_points) {
        int comp = static_cast<int>(i);
        if (comp == t) return true;
    }
    return false;
}
