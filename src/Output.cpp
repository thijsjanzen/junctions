//
//  Output.cpp
//  
//
//  Created by Thijs Janzen on 07/11/2017.
//
//

#include "Output.h"


void Output::update_inf(const std::vector<Fish_inf>& Pop) {
    double averageNumJunctions = 0;
    for(std::vector<Fish_inf>::const_iterator i = Pop.begin(); i != Pop.end(); ++i)
    {
        int numJ  = (int)(*i).chromosome1.size() - 2; //exclude the ends
        numJ += (int)(*i).chromosome2.size() - 2; //exclude the ends
        averageNumJunctions += numJ;
    }

    averageNumJunctions = 1.0 * averageNumJunctions / (2 * (int)Pop.size());

    avgJunctions.push_back(averageNumJunctions);

    return;
}

int countJunctions(const std::vector<bool>& B) {
    int numJunctions = 0;
    for(std::size_t i = 1; i < B.size(); ++i) {
        if(B[i] != B[i-1]) {
            numJunctions++;
        }
    }
    return numJunctions;
}


void Output::update_fin(const std::vector<Fish_fin>& Pop) {
    double averageNumJunctions = 0;
    for(auto i = Pop.begin(); i != Pop.end(); ++i)
    {
        int numJ = countJunctions((*i).chromosome1);
        averageNumJunctions += numJ;
        numJ = countJunctions((*i).chromosome2);
        averageNumJunctions += numJ;
    }

    averageNumJunctions = 1.0 * averageNumJunctions / (2*Pop.size());

    avgJunctions.push_back(averageNumJunctions);
    
    return;
}


std::vector<bool> detectJunctions(const std::vector<junction>& G,
                                  const std::vector<double>& markers) {
    std::vector<bool> output(markers.size());

    int j = 0;
    for(int i = 0; i < markers.size(); ++i) {
        double focalPos = markers[i];
        for(; j <= (G.size()-1); ++j) {
            double left = G[j].pos;
            double right = G[j+1].pos;
            if(left <= focalPos && right >= focalPos) {

                output[i] = ((bool)G[j].right);
                break;
            }
        }
    }
    return output;
}

void Output::detectNumJunctions(const std::vector<Fish_inf> &Pop,
                                const std::vector<double> &markers) {
    double averageNumJunctions = 0;
    for(auto i = Pop.begin(); i != Pop.end(); ++i) {
        std::vector<bool> genome1 = detectJunctions((*i).chromosome1, markers);
        averageNumJunctions += countJunctions(genome1);

        std::vector<bool> genome2 = detectJunctions((*i).chromosome2, markers);
        averageNumJunctions += countJunctions(genome2);
    }

    averageNumJunctions = 1.0 * averageNumJunctions / (2 * Pop.size()); //diploid
    avg_detected_Junctions.push_back(averageNumJunctions);
    return;
}
