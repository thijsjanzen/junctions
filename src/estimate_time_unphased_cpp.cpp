
#include <unistd.h>
#include <vector>
#include <cassert>
#include <complex.h>
#include <cmath>
#include <iostream>

#include "estimate_time_unphased_cpp.h"

#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

void force_output() {
  // Rcout << s << "\n";
  R_FlushConsole();
  R_ProcessEvents();
  R_CheckUserInterrupt();
  ::sleep(1);
}


struct chromosome {
  std::vector< int > states;
  std::vector< double > distances;

  chromosome(const NumericMatrix& anc_matrix,
             const NumericVector& loc) {
    states = std::vector<int>(loc.size(), 2);
    distances = std::vector<double>(loc.size(), 0.0);
    //Rcout << states.size() << " " << distances.size() << "\n"; force_output();

    for(int i = 0; i < distances.size(); ++i) {
      // assert(i < distance.size());
      // assert(i < anc_matrix.nrow());
      if (i > 0) distances[i] = loc[i] - loc[i - 1];
      if (anc_matrix(i, 0) == anc_matrix(i, 1)) {
        states[i] = anc_matrix(i, 0); // anc_matrix(i, 0) = [0, 1]
      }
      // Rcout << i << " " << states[i] << " " << distances[i] << "\n"; force_output();
    }
  }

  double calculate_likelihood(double t, int pop_size, double freq_ancestor_1);
};


// [[Rcpp::export]]
double estimate_time_unphased_cpp(const NumericMatrix& local_anc_matrix,
                                  const NumericVector& locations,
                                  int pop_size,
                                  double freq_ancestor_1,
                                  int lower_lim,
                                  int upper_lim,
                                  bool verbose) {

  chromosome focal_chrom(local_anc_matrix,
                         locations);
  return 0.0;
}

//' function to calculate log likelihood using cpp
//' @param local_anc_matrix local ancestry matrix
//' @param locations locations of markers
//' @param pop_size population size
//' @param freq_ancestor_1 frequency of the most common ancestor
//' @param t time
//' @export
// [[Rcpp::export]]
double loglikelihood_unphased_cpp(const NumericMatrix& local_anc_matrix,
                                  const NumericVector& locations,
                                  int pop_size,
                                  double freq_ancestor_1,
                                  double t) {

  Rcout << "loading chromosome\n"; //force_output();

  chromosome focal_chrom(local_anc_matrix,
                         locations);

  Rcout << "starting likelihood calculation\n"; force_output();

  double ll = focal_chrom.calculate_likelihood(t, pop_size, freq_ancestor_1);
  return ll;
}



//' function to calculate 7 states
//' @param t time
//' @param N pop size
//' @param d distance
//' @export
// [[Rcpp::export]]
std::vector< double > single_state_cpp(int t, int N, double d) {

  // I verified this with the synonymous R code, and it generates
  // the correct answer.
  // the cpp version is about 6 times faster:
  //                                            expr    min      lq     mean  median      uq     max neval cld
  // single_state(t = 100, N = 1000, d = 1e-05) 19.094 20.4415 22.72457 21.1695 21.9395 104.379   100   b
  // single_state_cpp(t = 100, N = 1000, d = 1e-05)  3.236  3.7970  4.32687  4.0560  4.3925  14.303   100  a

   double trans_matrix[7][7] =
    {{1.0 - 1.0 / (2*N) - 2 * d , 2 * d, 0, 0, 0, 1.0 / (2*N), 0},
    {1.0 / (2*N), 1 - 3 * 1.0 / (2*N) - d, d, 2 * 1.0 / (2*N), 0, 0, 0},
    {0, 2 * 1.0 / (2*N), 1 - 4 * 1.0 / (2*N), 0, 2 * 1.0 /(2*N), 0, 0},
    {0, 0, 0, 1 - 1.0 / (2*N) - d, d, 1.0 / (2*N), 0},
    {0, 0, 0, 2 * 1.0 / (2*N), 1 - 3 * 1.0 / (2*N), 0, 1.0 / (2*N)},
    {0, 0, 0, 0, 0, 1.0 - d, d},
    {0 ,0, 0, 0, 0, 1.0 / (2*N),  1 - 1.0 / (2*N)}};


  typedef SqMx<7, double> M7;
  M7 m(trans_matrix);
  auto m2 = m ^ t;
  std::vector< double > output = m2();

  return(output);
}



//' function to calculate prob
//' @param l left
//' @param r right
//' @param p freq
//' @param P seven states output
//' @export
// [[Rcpp::export]]
double get_prob_from_matrix_cpp(int left,
                                int right,
                                double p,
                                const std::vector<double>& P_) {
  // unphased probabilities!
  // left = marker on the left hand side [1 = PP, 2 = QQ, 3 = PQ or QP]
  // right = marker on the right hand side [1 = PP, 2 = QQ, 3 = PQ or QP]
  // p = frequency ancestor 1
  // P = states vector P_t^i

  left++; // to conform to R notation. // code below is a close copy of R code
  right++;

  std::vector<double> P(1, 0);
  P.insert(P.end(), P_.begin(), P_.end());

  // 0,0  0,1  0,2
  // 1,1  1,2  1,3

  double q  = 1 - p;
  double prob = 0;
  if (left == 1 && right == 1) {
    prob = (p * p) * (P[1] + P[4] + P[7]) +
      pow(p, 3) * (P[2] + P[5]) +
      pow(p, 4) * P[3] +
      p * P[6];
  }
  if (left == 1 && right == 2) {
    prob = p * q * (p * q * P[3] +
      (1.0 / 2) * P[5] +
      P[7]);
  }
  if (left == 1 && right == 3) {
    prob = p * q * (p * P[2] +
      2 * (p * p) * P[3] +
      (1.0 / 2) * P[4] +
      p * P[5]);
  }

  if (left == 2 && right == 1) {
    prob = p * q * (p * q * P[3] +
      (1.0 / 2) * P[5] +
      P[7]);
  }
  if (left == 2 && right == 2) {
    prob = (q * q) * (P[1] + P[4] + P[7]) +
      pow(q, 3) * (P[2] + P[5]) +
      (q * q) * P[3] +
      q * P[6];
  }
  if (left == 2 && right == 3) {
    prob = p * q * (q * P[2] +
      2 * (q * q) * P[3] +
      (1.0 / 2) * P[4] +
      q * P[5]);
  }

  if (left == 3 && right == 1) {
    prob = p * q * (p * P[2] +
      2 * (p * p) * P[3] +
      (1.0 / 2) * P[4] +
      p * P[5]);
  }
  if (left == 3 && right == 2) {
    prob = p * q * (q * P[2] +
      2 * (q * q) * P[3] +
      (1.0 / 2) * P[4] +
      q * P[5]);
  }
  if (left == 3 && right == 3) {
    prob = p * q * (2 * P[1] +
      P[2] +
      4 * p * q * P[3]);
  }
  return prob;
}



//' function to calculate ll
//' @param di distance
//' @param l left
//' @param r right
//' @param t time
//' @param pop_size population size
//' @param freq_ancestor_1 freq ancestor
//' @param condition conditioning
//' @export
// [[Rcpp::export]]
double calc_ll(double di,
               double l,
               double r,
               double t,
               int pop_size,
               double freq_ancestor_1,
               bool condition) {
  if (di < 0)
    return(-1e20);

  std::vector<  double > seven_states = single_state_cpp(t, pop_size, di);
  for(int i = 0; i < seven_states.size(); ++i) {
    Rcout << seven_states[i] << " ";
  }
  Rcout << "\n";

  std::vector<  double > probs(3);
  double sum_prob = 0.0;
  for(int i = 0; i < 3; ++i) {
    probs[i] = get_prob_from_matrix_cpp(l,
                                        i,
                                        freq_ancestor_1,
                                        seven_states);
    sum_prob += probs[i];
    Rcout << probs[i] << " ";
  }
  Rcout << "\n";

  double focal_prob = probs[r];
  Rcout << focal_prob << " " << sum_prob << "\n";

  if (condition == true)
    focal_prob *= 1.0 / sum_prob;

  return(log(focal_prob));
}

double chromosome::calculate_likelihood(double t,
                                        int pop_size,
                                        double freq_ancestor_1) {

  if (t < 1)
    return(-1e20);
  if (pop_size < 2)
    return(-1e20);
  if (freq_ancestor_1 >= 1)
    return(-1e20);
  if (freq_ancestor_1 <= 0)
    return(-1e20);


  double di = distances[1];
  double l = states[0];
  double r = states[1];

  //  double ll = calc_ll(di, l, r, t, pop_size, freq_ancestor_1, false);

  std::vector< double> ll(distances.size());
  ll[0] = calc_ll(di, l, r, t, pop_size, freq_ancestor_1, false);
  for(int i = 1; i < distances.size(); ++i) {
    //  Rcout << i << "\n"; force_output();
    di = distances[i];
    l = states[i - 1];
    r = states[i];
    ll[i] = calc_ll(di, l, r, t, pop_size, freq_ancestor_1, true);
    // std::cout << di << " " << l << " " << r << " " << ll[i] << "\n";
  }

  return(std::accumulate(ll.begin(), ll.end(), 0.0));
}

