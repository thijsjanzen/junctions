//
//  random_functions.cpp
//
//
//  Created by Thijs Janzen on 05/03/2018.
//
//

#include "random_functions.h"
#include <Rcpp.h>
#include <random>
using namespace Rcpp;

std::random_device rd;
std::mt19937 rndgen(rd());  //< The one and only random number generator

int random_number(int n)    {
    return std::uniform_int_distribution<> (0, n-1)(rndgen);
}

double uniform()    {
    return std::uniform_real_distribution<>(0, 1.0)(rndgen);
}

int poisson(double lambda)    {
    return std::poisson_distribution<int>(lambda)(rndgen);
}

void set_seed(unsigned seed)    {
    std::mt19937 new_randomizer(seed);
    rndgen = new_randomizer;
}