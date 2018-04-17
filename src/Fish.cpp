//
//  Fish.cpp
//  
//
//  Created by Thijs Janzen on 07/11/2017.
//
//

#include "Fish.h"
#include "random_functions.h"

double getRecomPos() {
    double pos = uniform();
    
    return pos;
}

int getRecomPos(int L) {
    int pos = -100;
    int index = random_number(L);
    // exclude the ends of the chromosome
    while (index == 0 || index == L) {
        index = random_number(L);
    }
    pos = index;

    return pos;
}

void Recombine_inf(std::vector<junction>& offspring,
               std::vector<junction> chromosome1,
               std::vector<junction> chromosome2,
               double MORGAN)  {

    int numRecombinations = poisson(MORGAN);

    if (numRecombinations == 0) {
        offspring.insert(offspring.end(),
                         chromosome1.begin(),
                         chromosome1.end());

        return;
    }

    //if the number of recombinations is larger than 1, we need some more complicated (slower) code:

    std::vector<double> recomPos(numRecombinations, 0);
    for(int i = 0; i < numRecombinations; ++i) {
        recomPos[i] = getRecomPos();
    }
    std::sort(recomPos.begin(), recomPos.end() );
    recomPos.erase(std::unique(recomPos.begin(), recomPos.end()), recomPos.end());

    while (recomPos.size() < numRecombinations) {
        double pos = getRecomPos();
        recomPos.push_back(pos);
        // sort them, in case they are not sorted yet
        // we need this to remove duplicates, and later
        // to apply crossover
        std::sort(recomPos.begin(), recomPos.end() );
        // remove duplicate recombination sites
        recomPos.erase(std::unique(recomPos.begin(), recomPos.end()), recomPos.end());
    }

    std::vector< junction > toAdd; //first create junctions on exactly the recombination positions
    for(int i = 0; i < recomPos.size(); ++i) {
        junction temp;
        temp.pos = recomPos[i];
        toAdd.push_back(temp);
    }

    //for(int i = 1; i < chromosome1.size(); ++i) {
    for(auto i = (chromosome1.begin()+1); i != chromosome1.end(); ++i) {
        double leftpos = (*(i-1)).pos;
        double rightpos = (*i).pos;

        for(int j = 0; j < recomPos.size(); ++j) {
            if(recomPos[j] > leftpos) {
                if(recomPos[j] < rightpos) {
                    if(j % 2 == 0) { //even, so chrom1 = L, chrom2 = R
                        toAdd[j].left = (*i).left;
                    } else { //uneven so chrom1 = R, chrom2 = L
                        toAdd[j].right = (*i).left;
                    }
                }
            }
        }
    }

    for(auto i = (chromosome2.begin()+1); i != chromosome2.end(); ++i) {
        double leftpos = (*(i-1)).pos;
        double rightpos = (*i).pos;

        for(int j = 0; j < recomPos.size(); ++j) {
            if(recomPos[j] > leftpos) {
                if(recomPos[j] < rightpos) {
                    if(j % 2 == 0) { //even, so chrom1 = L, chrom2 = R
                        toAdd[j].right = (*i).left;
                    } else { //uneven so chrom1 = R, chrom2 = L
                        toAdd[j].left = (*i).left;
                    }
                }
            }
        }
    }

    for(int i = 0; i < toAdd.size(); ++i) {
        if(toAdd[i].left != toAdd[i].right) {
            offspring.push_back(toAdd[i]);
        }
    }

    //now we have to add the other junctions from chrom1 and chrom2.
    double leftpos = 0;
    double rightpos = 0;


    for(int i = 0; i < (recomPos.size() + 1); ++i) {
        rightpos = 1.0;
        if(i < recomPos.size()) rightpos = recomPos[i];

        if(i % 2 == 0) { //even, so take from chromosome 1
            for(std::vector<junction>::iterator it = chromosome1.begin(); it != chromosome1.end(); ++it) {
                if((*it).pos >= leftpos && (*it).pos <= rightpos) {
                    offspring.push_back((*it));
                }
                if((*it).pos > rightpos) { //we are past the recombination section
                    break;
                }
            }
        }

        if(i % 2 == 1) { //odd, so take from chromosome 2
            for(std::vector<junction>::iterator it = chromosome2.begin(); it != chromosome2.end(); ++it) {
                if((*it).pos >= leftpos && (*it).pos <= rightpos) {
                    offspring.push_back((*it));
                }
                if((*it).pos > rightpos) { //we are past the recombination section
                    break;
                }
            }
        }

        //move forward
        leftpos = rightpos;
    }

    std::sort(offspring.begin(), offspring.end());
    offspring.erase(std::unique(offspring.begin(), offspring.end()), offspring.end());
    
    return;
}

Fish_inf mate_inf(const Fish_inf& A, const Fish_inf& B, double numRecombinations)
{
    Fish_inf offspring;
    offspring.chromosome1.clear();
    offspring.chromosome2.clear(); //just to be sure.

    //first the father chromosome
    int event = random_number(2);
    switch(event) {
        case 0:  {
            Recombine_inf(offspring.chromosome1, A.chromosome1, A.chromosome2, numRecombinations);
            break;
        }
        case 1: {
            Recombine_inf(offspring.chromosome1, A.chromosome2, A.chromosome1, numRecombinations);
            break;
        }
    }


    //then the mother chromosome
    event = random_number(2);
    switch(event) {
        case 0:  {
            Recombine_inf(offspring.chromosome2, B.chromosome1, B.chromosome2, numRecombinations);
            break;
        }
        case 1: {
            Recombine_inf(offspring.chromosome2, B.chromosome2, B.chromosome1, numRecombinations);
            break;
        }
    }

    return offspring;
}

void Recombine_fin(std::vector<bool>* offspring,
               std::vector<bool> chromosome1,
               std::vector<bool> chromosome2,
               double numberRecombinations)  {
    // we have a decimal value, the fractional part is interpreted as
    // a probability of one extra recombination, note that we do NOT
    // assume a poisson distributed number of recombinations
    // this has been proven inaccurate, and tends to overestimate
    // the number of meioses without recombination

    numberRecombinations = poisson(numberRecombinations);

    // if there are not recombinations, preliminary exit
    if (numberRecombinations == 0) {
        //offspring = &chromosome1;
        offspring->insert(offspring->end(),
                          chromosome1.begin(),
                          chromosome1.end() );
        return;
    }

    std::vector<int> recomPos;
    // store L, so we avoid repeated calls of the function .size()
    int L = static_cast<int>(chromosome1.size());

    while (recomPos.size() < numberRecombinations) {
        int pos = getRecomPos(L);
        recomPos.push_back(pos);
        // sort them, in case they are not sorted yet
        // we need this to remove duplicates, and later
        // to apply crossover
        std::sort(recomPos.begin(), recomPos.end());
        // remove duplicate recombination sites
        std::vector<int>::iterator last = std::unique(recomPos.begin(), recomPos.end());
        recomPos.erase(last, recomPos.end());
    }

    // used to track which chromosome was used
    // during the last recombination event
    int order = 0;
    int start = 0;

    for (std::size_t i = 0; i < recomPos.size(); ++i) {
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
              double numberRecombinations) {

    Fish_fin offspring;
    offspring.chromosome1.clear();
    offspring.chromosome2.clear();  // just to be sure.

    // random order or in other words,
    // we randomly select 1 of 2 produced chromosomes during recombination
    if (uniform() < 0.5) {
        Recombine_fin(&offspring.chromosome1,
                  A.chromosome1, A.chromosome2,
                  numberRecombinations);
    } else {
        Recombine_fin(&offspring.chromosome1,
                  A.chromosome2, A.chromosome1,
                  numberRecombinations);
    }

    if (uniform() < 0.5) {
        Recombine_fin(&offspring.chromosome2,
                  B.chromosome1, B.chromosome2,
                  numberRecombinations);
    } else {
        Recombine_fin(&offspring.chromosome2,
                  B.chromosome2, B.chromosome1,
                  numberRecombinations);
    }
    return offspring;
}
