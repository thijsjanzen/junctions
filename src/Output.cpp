//
//  Output.cpp
//
//
//  Created by Thijs Janzen on 07/11/2017.
//
//

#include "Output.h"
#include "Rcpp.h"
using namespace Rcpp;

void force_output() {
    // Rcout << s << "\n";
    static std::chrono::milliseconds timespan(100);
    std::this_thread::sleep_for(timespan);
    R_FlushConsole();
    R_ProcessEvents();
    R_CheckUserInterrupt();
}

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
    for(unsigned int i = 1; i < B.size(); ++i) {
        if(B[i] != B[i-1]) {
            numJunctions++;
        }
    }
    return numJunctions;
}


void Output::update_fin(const std::vector<Fish_fin>& Pop) {
    double averageNumJunctions = 0;
    for(auto i = Pop.begin(); i != Pop.end(); ++i) {
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

    unsigned int j = 0;
    for(unsigned int i = 0; i < markers.size(); ++i) {
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

std::vector<int> detect_ancestry(const std::vector< junction >& G,
                                 const std::vector< double >& markers) {
    std::vector<int> output(markers.size());
   // Rcout << "detect_ancestry\n"; force_output();
    int j = 0;
    for(int i = 0; i < static_cast<int>(markers.size()); ++i) {
        double focalPos = markers[i];
        for(; j <= (G.size()-1); ++j) {
            double left = G[j].pos;
            double right = G[j+1].pos;
            if(left <= focalPos && right >= focalPos) {
                output[i] = ((int)G[j].right);
                break;
            }
        }
        j-=5; //just in case
        if(j < 0) j = 0;
    }
    return output;
}

void Output::update_unphased(const std::vector< Fish_inf >& Pop,
                    size_t t,
                    bool record_true_junctions,
                    double morgan,
                    size_t num_indiv) {

    for(unsigned int i = 0; i < num_indiv; ++i) {
        std::vector<int> chrom1 = detect_ancestry(Pop[i].chromosome1, markers);
        std::vector<int> chrom2 = detect_ancestry(Pop[i].chromosome2, markers);

        for(unsigned int j = 0; j < markers.size(); ++j) {
            std::vector<double> to_add(5); // = {t, i, markers[j], chrom1[j], chrom2[j]};
            to_add[0] = t;
            to_add[1] = i; //individual
            to_add[2] = markers[j] * morgan;
            to_add[3] = chrom1[j];
            to_add[4] = chrom2[j];
            results.push_back(to_add);
        }

        if(record_true_junctions) {
            int true_junct_chrom_1  = (int)Pop[i].chromosome1.size() - 2; //exclude the ends
            int true_junct_chrom_2  = (int)Pop[i].chromosome2.size() - 2; //exclude the ends
            std::vector< double > to_add_true(4);
            to_add_true[0] = t;
            to_add_true[1] = i;
            to_add_true[2] = true_junct_chrom_1;
            to_add_true[3] = true_junct_chrom_2;
            true_results.push_back(to_add_true);
        }
    }
    return;
}


int detect_junctions(const Fish_inf& indiv,
                     const std::vector<double> &markers,
                     double& avg_heterozygosity) {

    // we need to find if at specific markers, they are homozygous or heterozygous
    std::vector<bool> chrom1 = detectJunctions(indiv.chromosome1, markers);
    std::vector<bool> chrom2 = detectJunctions(indiv.chromosome2, markers);


    std::vector<int> genotypes(chrom1.size(), -1);
    for(unsigned int i = 0; i < chrom1.size(); ++i) {
        // 0 = homozygous parent 0
        // 1 = heterozygous
        // 2 = homozygous parent 1

        if (chrom1[i] != chrom2[i]) { // heterozygous, e.g. 0/1 or 1/0
            genotypes[i] = 1;
        } else {
            if (chrom1[i] == 0) {  // homozygous 0 : 0/0
                genotypes[i] = 0;
            } else {               // homozygous 1 : 1/1
                genotypes[i] = 2;
            }
        }

       /* if(chrom1[i] == 0 && chrom2[i] == 0) {
            genotypes[i] = 0;
        }
        if(chrom1[i] == 0 && chrom2[i] == 1) {
            genotypes[i] = 1;
        }
        if(chrom1[i] == 1 && chrom2[i] == 0) {
            genotypes[i] = 1;
        }
        if(chrom1[i] == 1 && chrom2[i] == 1) {
            genotypes[i] = 2;
        }*/
    }

    int number_of_junctions = 0;
    int number_heterozygous = 0;

    // zero entry is not included in the loop
    if (genotypes[0] == 1) number_heterozygous++;

    for (unsigned int i = 1; i < genotypes.size(); ++i) {
        if(genotypes[i] != -1 && genotypes[i-1] != -1) {
            if(genotypes[i] != genotypes[i-1]) number_of_junctions++;
        }
        if(genotypes[i] == 1) number_heterozygous++;
    }
    avg_heterozygosity += 1.0 * number_heterozygous / markers.size();

    return number_of_junctions;
}


void Output::detect_junctions_backcross(const std::vector< Fish_inf > &Pop,
                                        const std::vector<double> &markers) {

    double average_detected_junctions = 0;

    std::vector<int> J;
    double avg_heterozygosity = 0.0;
    for(auto i = Pop.begin(); i != Pop.end(); ++i) {
        int dJ = detect_junctions((*i), markers, avg_heterozygosity);
        average_detected_junctions += dJ;
        J.push_back(dJ);
    }
    junction_dist.push_back(J);

    average_detected_junctions = 1.0 * average_detected_junctions / (2 * Pop.size()); //diploid
    avg_detected_Junctions.push_back(average_detected_junctions);
    avg_hetero.push_back(avg_heterozygosity / Pop.size());
    return;
}