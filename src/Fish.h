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
    long double pos;
    int left, right;

    junction()  {
    }

    junction(long double loc, int A, int B) : pos(loc), left(A), right(B) {
    }

    junction(const junction& other) {
        pos = other.pos;
        left = other.left;
        right = other.right;
    }

    bool operator ==(const junction& other) const {
        if(pos != other.pos) return false;
        if(left != other.left) return false;
        if(right != other.right) return false;

        return true;
    }

    bool operator <(const junction& other) const {
        return(pos < other.pos);
    }

    bool operator !=(const junction& other) const {
        return( !( (*this) == other) );
    }
};


struct Fish_inf {
    std::vector< junction > chromosome1;
    std::vector< junction > chromosome2;

    Fish_inf()
    {}

    Fish_inf(int initLoc)    {
        junction left = junction(0.0, -1, initLoc);
        junction right = junction(1, initLoc, -1);
        chromosome1.push_back( left  );
        chromosome1.push_back( right );
        chromosome2.push_back( left  );
        chromosome2.push_back( right );
    }

    Fish_inf(const std::vector<junction>& A,
             const std::vector<junction>& B)    {
        chromosome1 = A;
        chromosome2 = B;
    }
};

struct Fish_fin  {
    std::vector<bool> chromosome1;
    std::vector<bool> chromosome2;

    Fish_fin() {
    }

    // constructor that sets all genome elements to "initLoc"
    Fish_fin(const bool initLoc, const int genomeSize) {
        for ( int i = 0; i < genomeSize; ++i ) {
            chromosome1.push_back(initLoc);
            chromosome2.push_back(initLoc);
        }
    }

    // copy constructor
    Fish_fin(const std::vector<bool>& A, const std::vector<bool>& B) {
        chromosome1 = A;
        chromosome2 = B;
    }
};


Fish_fin mate_fin(const Fish_fin& A, const Fish_fin& B,
          double numRecombinations);

Fish_inf mate_inf(const Fish_inf& A, const Fish_inf& B,
                  double numRecombinations);

long double getRecomPos();
int getRecomPos(int L);




#endif /* Fish_hpp */
