//
//  Fish.cpp
//
//
//  Created by Thijs Janzen on 07/11/2017.
//
//

#include "Fish.h"
#include "random_functions.h"

#include <Rcpp.h>
using namespace Rcpp;

void add(std::vector< junction>& offspring,
         const junction& new_junction) {

    if (offspring.empty()) {
        offspring.push_back(new_junction);
        return;
    }


    if (new_junction.pos > offspring.back().pos &&
        new_junction.right != offspring.back().right) {
        offspring.push_back(new_junction);
    }
    return;
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
void do_recombination(std::vector<junction>& offspring,
                      const std::vector<junction>& chromosome1,
                      const std::vector<junction>& chromosome2,
                      std::vector<double>& recomPos) {

    /*
     chrom 1:  [0    1
     0.5  0
     1.0 -1]

     chrom 2:  [0    0
     0.5  1
     1.0 -1]

     recomPos  = {0.75, 0.75};

     new_chrom = [0      1
     0.5    0
     1.0   -1]

     // std::lower_bound()

     */

    std::vector < std::vector<junction>::const_iterator > iters =
        { chromosome1.begin(), chromosome2.begin() };

 //   recomPos.push_back(1.0); // for completeness

    int index = 0;
    int recompos_cnt = 0;

    while(true) {

        if ( iters[index]->pos > recomPos[recompos_cnt]  ) {
            // encountered junction point
            // create junction
            index = 1 - index;
            while( iters[index]->pos < recomPos[recompos_cnt]) {
                iters[index]++;
            }

            auto prev_iter = iters[index];
            prev_iter--;
            junction new_junction(recomPos[recompos_cnt], prev_iter->right);
            add(offspring, new_junction);

            recompos_cnt++;
        } else {
            add(offspring, (*iters[index]));
            iters[index]++;
        }

        if (offspring.back().right == -1) {
            break;
        }
    }

    return;
}

std::vector<junction> recombine_new(const std::vector<junction>& chromosome1,
                                    const std::vector<junction>& chromosome2,
                                    const std::vector<double>& recom_positions) {

    static thread_local auto tl_go = decltype(chromosome1){};
    assert(!chromosome1.empty());    // not strictly enforced by code
    assert(!chromosome2.empty());    // not strictly enforced bu code

    // we need something that is cheaply swappable:
    auto* g1 = &chromosome1;
    auto* g2 = &chromosome2;
    auto& go = tl_go;   // offspring genome: recycle what's already there...
    go.clear();

    // predicate for lower_bound
    auto less = [](const auto& j, double p) { return j.pos < p; };

    // helper lambda to get the value just *before* it.
    // we store the value to the right of a recombination-point but we need the value to the left:
    auto value_at = [](auto begin, auto it) { return (begin != it) ? (it - 1)->right : -1; };

    double left_pos = 0.0;
    auto go_val = -1;
    for (auto right_pos : recom_positions) {
        auto it = std::lower_bound(g1->cbegin(), g1->cend(), left_pos, less);
        auto last = std::lower_bound(it, g1->cend(), right_pos, less);
        // [g1.first, it) : part of the genome *before* left_pos.
        // [it, last) : part of the genome *after or equal to* left_pos but *before* right_pos.
        auto g1_val = value_at(g1->cbegin(), it);
        if (g1_val != go_val) {
            if (it == last || it->pos != left_pos) {
                go.emplace_back(left_pos, g1_val);   // insert change to match
            }
            else {
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

std::vector<double> generate_recomPos(int number_of_recombinations) {

    std::vector<double> recomPos(number_of_recombinations, 0);
    for(int i = 0; i < number_of_recombinations; ++i) {
        recomPos[i] = uniform();
    }
    std::sort(recomPos.begin(), recomPos.end());
    if (recomPos.size() != number_of_recombinations) {
        stop("mismatch\n");
    }
    recomPos.push_back(1.0);

    return recomPos;
}

void Recombine_inf(      std::vector<junction>& offspring,
                     const std::vector<junction>& chromosome1,
                     const std::vector<junction>& chromosome2,
                     double MORGAN)  {

    int numRecombinations = poisson(MORGAN);
   // Rcout << numRecombinations << "\n";

    if (numRecombinations == 0) {
        offspring.insert(offspring.end(),
                         chromosome1.begin(),
                         chromosome1.end());

        return;
    }

    std::vector<double> recomPos = generate_recomPos(numRecombinations);

    offspring = recombine_new(chromosome1,
                              chromosome2,
                              recomPos);
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

    for (int i = 0; i < recomPos.size(); ++i) {
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
/////////////////////////////////////////////
////////////// Member functions
/////////////////////////////////////////////


junction::junction() {
}

junction::junction(double loc, int A) : pos(loc), right(A) {
}

junction::junction(const junction& other) {
    pos = other.pos;
    right = other.right;
}

bool junction::operator <(const junction& other) const {
    return(pos < other.pos);
}

Fish_inf::Fish_inf(){

}

Fish_inf::Fish_inf(int initLoc)    {
    junction left = junction(0.0, initLoc);
    junction right = junction(1, -1);
    chromosome1.push_back( left  );
    chromosome1.push_back( right );
    chromosome2.push_back( left  );
    chromosome2.push_back( right );
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

bool is_in_time_points(int t,
                       const Rcpp::NumericVector & time_points) {
    for (auto it = time_points.begin(); it != time_points.end(); ++it) {
        if((*it) == t) return true;
    }
    return false;
}
