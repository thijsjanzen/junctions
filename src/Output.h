//
//  Output.hpp
//  
//
//  Created by Thijs Janzen on 07/11/2017.
//
//

#ifndef Output_hpp
#define Output_hpp

#include <stdio.h>
#include "Fish.h"
#include <vector>

struct Output {
    std::vector<double> avgJunctions;
    std::vector<double> avg_detected_Junctions;
    void update_inf(const std::vector< Fish_inf >& Pop);
    void update_fin(const std::vector< Fish_fin >& Pop);

    void detectNumJunctions(const std::vector<Fish_inf> &Pop,
                            const std::vector<double> &markers);
};

#endif /* Output_hpp */
