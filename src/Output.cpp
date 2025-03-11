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
//  Output.cpp
//
//
//  Created by Thijs Janzen on 07/11/2017.
//
//


#include "Output.h" // NOLINT [build/include_subdir]
#include "Rcpp.h"
#include <vector>

void Output::update_inf(const std::vector<Fish_inf>& Pop) {
    double averageNumJunctions = 0;
    for (const auto& i : Pop) {
        //  exclude the ends:
        int numJ  = static_cast<int>(i.chromosome1.size()) - 2 +
                    static_cast<int>(i.chromosome2.size()) - 2;
        averageNumJunctions += numJ;
    }
    averageNumJunctions = 1.0 * averageNumJunctions / (2 * Pop.size());
    avgJunctions.push_back(averageNumJunctions);
    return;
}

int countJunctions(const std::vector<bool>& B) {
    int numJunctions = 0;
    for (unsigned int i = 1; i < B.size(); ++i) {
        if (B[i] != B[i-1]) {
            numJunctions++;
        }
    }
    return numJunctions;
}


void Output::update_fin(const std::vector<Fish_fin>& Pop) {
    double averageNumJunctions = 0;
    for (const auto& i : Pop) {
        averageNumJunctions += countJunctions(i.chromosome1);;
        averageNumJunctions += countJunctions(i.chromosome2);
    }

    averageNumJunctions = 1.0 * averageNumJunctions / (2 * Pop.size());
    avgJunctions.push_back(averageNumJunctions);
    return;
}


std::vector<bool> detectJunctions(const std::vector<junction>& G,
                                  const std::vector<double>& markers) {
    std::vector<bool> output(markers.size());

    size_t j = 0;
    for (size_t i = 0; i < markers.size(); ++i) {
        double focalPos = markers[i];
        for (; j <= (G.size() - 1); ++j) {
            double left = G[j].pos;
            double right = G[j + 1].pos;
            if (left <= focalPos && right >= focalPos) {
                output[i] = static_cast<bool>(G[j].right);
                break;
            }
        }
    }
    return output;
}

void Output::detectNumJunctions(const std::vector<Fish_inf> &Pop,
                                const std::vector<double> &markers) {
    double averageNumJunctions = 0;
    for (const auto& i : Pop) {
        std::vector<bool> genome1 = detectJunctions(i.chromosome1, markers);
        averageNumJunctions += countJunctions(genome1);

        std::vector<bool> genome2 = detectJunctions(i.chromosome2, markers);
        averageNumJunctions += countJunctions(genome2);
    }

    averageNumJunctions = 1.0 * averageNumJunctions / (2 * Pop.size());
    avg_detected_Junctions.push_back(averageNumJunctions);
    return;
}

std::vector<int> detect_ancestry(const std::vector< junction >& G,
                                 const std::vector< double >& markers) {
    std::vector<int> output(markers.size());
    int j = 0;
    for (int i = 0; i < static_cast<int>(markers.size()); ++i) {
        double focalPos = markers[i];
        for (; j <= (G.size() - 1); ++j) {
            double left = G[j].pos;
            double right = G[j + 1].pos;
            if (left <= focalPos && right >= focalPos) {
                output[i] = static_cast<int>(G[j].right);
                break;
            }
        }
        j -= 5;   // just in case
        if (j < 0) j = 0;
    }
    return output;
}

void Output::update_unphased(const std::vector< Fish_inf >& Pop,
                             size_t t,
                             bool record_true_junctions,
                             double morgan,
                             size_t num_indiv) {
    for (size_t i = 0; i < num_indiv; ++i) {
        std::vector<int> chrom1 = detect_ancestry(Pop[i].chromosome1, markers);
        std::vector<int> chrom2 = detect_ancestry(Pop[i].chromosome2, markers);
        for (size_t j = 0; j < markers.size(); ++j) {
            std::vector<double> to_add(5);
            to_add[0] = t;
            to_add[1] = i;   // individual
            to_add[2] = markers[j] * morgan;
            to_add[3] = chrom1[j];
            to_add[4] = chrom2[j];
            results.push_back(to_add);
        }

        if (record_true_junctions) {
            std::vector< double > to_add_true(4);
            to_add_true[0] = t;
            to_add_true[1] = i;
            to_add_true[2] = Pop[i].chromosome1.size() - 2;
            to_add_true[3] = Pop[i].chromosome2.size() - 2;
            true_results.push_back(to_add_true);
        }
    }
    return;
}

int detect_junctions(const Fish_inf& indiv,
                     const std::vector<double> &markers,
                     double* avg_heterozygosity) {
    // we need to find if at specific markers,
    // they are homozygous or heterozygous
    std::vector<bool> chrom1 = detectJunctions(indiv.chromosome1, markers);
    std::vector<bool> chrom2 = detectJunctions(indiv.chromosome2, markers);

    std::vector<int> genotypes(chrom1.size(), -1);
    for (unsigned int i = 0; i < chrom1.size(); ++i) {
        // 0 = homozygous parent 0
        // 1 = heterozygous
        // 2 = homozygous parent 1

        if (chrom1[i] != chrom2[i]) {   //  heterozygous, e.g. 0/1 or 1/0
            genotypes[i] = 1;
        } else {
            if (chrom1[i] == 0) {  // homozygous 0 : 0/0
                genotypes[i] = 0;
            } else {               // homozygous 1 : 1/1
                genotypes[i] = 2;
            }
        }
    }

    int number_of_junctions = 0;
    int number_heterozygous = 0;

    // zero entry is not included in the loop
    if (genotypes[0] == 1) number_heterozygous++;

    for (size_t i = 1; i < genotypes.size(); ++i) {
        if (genotypes[i] != -1 && genotypes[i-1] != -1) {
            if (genotypes[i] != genotypes[i-1]) number_of_junctions++;
        }
        if (genotypes[i] == 1) number_heterozygous++;
    }
    *avg_heterozygosity += 1.0 * number_heterozygous / markers.size();

    return number_of_junctions;
}


void Output::detect_junctions_backcross(const std::vector< Fish_inf > &Pop,
                                        const std::vector<double> &markers) {
    double average_detected_junctions = 0;

    std::vector<int> J;
    double avg_heterozygosity = 0.0;
    for (const auto& i : Pop) {
        int dJ = detect_junctions(i, markers, &avg_heterozygosity);
        average_detected_junctions += dJ;
        J.push_back(dJ);
    }
    junction_dist.push_back(J);

    average_detected_junctions =
             1.0 * average_detected_junctions / (2 * Pop.size());   //  diploid
    avg_detected_Junctions.push_back(average_detected_junctions);
    avg_hetero.push_back(avg_heterozygosity / Pop.size());
    return;
}
