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

struct junction {
    double pos;
    int left, right;

    junction();
    junction(double loc, int A, int B);
    junction(const junction& other);

    bool operator <(const junction& other) const;
};


struct Fish_inf {
    std::vector< junction > chromosome1;
    std::vector< junction > chromosome2;

    Fish_inf();

    Fish_inf(int initLoc);
};

struct Fish_fin  {
    std::vector<bool> chromosome1;
    std::vector<bool> chromosome2;

    Fish_fin();

    Fish_fin(const bool initLoc, const int genomeSize);
};


Fish_fin mate_fin(const Fish_fin& A, const Fish_fin& B,
          double numRecombinations);

Fish_inf mate_inf(const Fish_inf& A, const Fish_inf& B,
                  double numRecombinations);

long double getRecomPos();
int getRecomPos(int L);
#endif /* Fish_hpp */
