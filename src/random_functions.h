//
//  random_functions.hpp
//
//
//  Created by Thijs Janzen on 05/03/2018.
//
//

#ifndef random_functions_hpp
#define random_functions_hpp

#include <random>
#include <vector>

double uniform();
int random_number(int n);
int poisson(double lambda);
void set_seed(unsigned seed);

std::vector<double> generate_random_markers(int number_of_markers);

#endif /* random_functions_hpp */
