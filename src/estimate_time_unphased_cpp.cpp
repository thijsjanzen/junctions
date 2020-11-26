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

struct chromosome {
  std::vector< int > states;
  std::vector< double > distances;

  chromosome(const NumericMatrix& anc_matrix,
             const NumericVector& loc) {
    states = std::vector<int>(3, loc.size());
    distances = std::vector<double>(0.0, loc.size());

    for(int i = 0; i < loc.size(); ++i) {
      if (i > 0) distances[i] = loc[i] - loc[i - 1];
      if (anc_matrix(i, 0) == anc_matrix(i, 1) && anc_matrix(i, 0) == 0) {
        states[i] = 1 + anc_matrix(i, 0); // anc_matrix(i, 0) = [0, 1], labels are [1, 2]
      }
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

  chromosome focal_chrom(local_anc_matrix,
                         locations);

  double ll = focal_chrom.calculate_likelihood(t, pop_size, freq_ancestor_1);
  return ll;
}



void mul(std::vector< std::vector< double > >& a,
         const std::vector< std::vector< double > >& b)
{
  std::vector< double > filler(a[0].size());
  std::vector< std::vector< double > > res(a.size(), filler);

  for (int i = 0; i < a.size(); i++) {
    for (int j = 0; j < a.size(); j++) {
      for (int k = 0; k < a.size(); k++)
      {
        res[i][j] += a[i][k] * b[k][j];
      }
    }
  }

  a = res;
  return;
}

void matrix_pow(std::vector< std::vector< double > >& trans_matrix,
                int t) {


  std::vector< double > filler(0, trans_matrix[0].size());
  std::vector< std::vector< double > > res(trans_matrix[0].size(), filler);
  // diagonal matrix
  for(int i = 0; i < trans_matrix[0].size(); ++i) {
    for(int j = 0; j < trans_matrix[0].size(); ++j) {
      res[i][j] = (i == j);
    }
  }

  while (t > 0) {
    if (t % 2 == 0)
    {
      mul(trans_matrix, trans_matrix);
      t /= 2;
    }
    else {
      mul(res, trans_matrix);
      t--;
    }
  }
}


std::vector< double > calc_mult_matrix(const std::vector< std::vector< double > >& trans_matrix,
                        int t) {

  std::vector< std::vector< double > > out_matrix = trans_matrix;
  matrix_pow(out_matrix, t);
  return(out_matrix[0]);  // same as multiplying with {1, 0, 0, 0, 0, 0, 0}
}

std::vector< double > single_state(int t, int N, double d) {
  std::vector<double> filler(7);
  std::vector< std::vector< double > > trans_matrix(7, filler);
  trans_matrix[0] = {1.0 - 1.0 / (2*N) - 2 * d , 2 * d, 0, 0, 0, 1.0 / (2*N), 0};
  trans_matrix[1] = {1.0 / (2*N), 1 - 3 * 1.0 / (2*N) - d, d, 2 * 1.0 / (2*N), 0, 0, 0};
  trans_matrix[2] = {0, 2 * 1.0 / (2*N), 1 - 4 * 1.0 / (2*N), 0, 2 * 1.0 /(2*N), 0, 0};
  trans_matrix[3] = {0, 0, 0, 1 - 1.0 / (2*N) - d, d, 1.0 / (2*N), 0};
  trans_matrix[4] = {0, 0, 0, 2 * 1.0 / (2*N), 1 - 3 * 1.0 / (2*N), 0, 1.0 / (2*N)};
  trans_matrix[5] = {0, 0, 0, 0, 0, 1.0 - d, d};
  trans_matrix[6] = {0 ,0, 0, 0, 0, 1.0 / (2*N),  1 - 1.0 / (2*N)};

  std::vector<double> result = calc_mult_matrix(trans_matrix, t);
  return(result);
}

double get_prob_from_matrix(int left,
                            int right,
                            double p,
                            const std::vector<double>& P) {
// unphased probabilities!
// left = marker on the left hand side [1 = PP, 2 = QQ, 3 = PQ or QP]
// right = marker on the right hand side [1 = PP, 2 = QQ, 3 = PQ or QP]
// p = frequency ancestor 1
// P = states vector P_t^i

  double q  = 1 - p;
  double prob = 0;
  if (left == 1 && right == 1) {
    prob = (p * p) * (P[1] + P[4] + P[7]) +
      (p * p * p) * (P[2] + P[5]) +
      (p * p * p * p) * P[3] +
      p * P[6];
  }
  if (left == 1 && right == 2) {
    prob = p * q * (p * q * P[3] +
      (1 / 2) * P[5] +
      P[7]);
  }
  if (left == 1 && right == 3) {
    prob = p * q * (p * P[2] +
      2 * (p * p) * P[3] +
      (1 / 2) * P[4] +
      p * P[5]);
  }

  if (left == 2 && right == 1) {
    prob = p * q * (p * q * P[3] +
      (1 / 2) * P[5] +
      P[7]);
  }
  if (left == 2 && right == 2) {
    prob = (q * q) * (P[1] + P[4] + P[7]) +
      (q * q * q) * (P[2] + P[5]) +
      (q * q) * P[3] +
      q * P[6];
  }
  if (left == 2 && right == 3) {
    prob = p * q * (q * P[2] +
      2 * (q * q) * P[3] +
      (1 / 2) * P[4] +
      q * P[5]);
  }

  if (left == 3 && right == 1) {
    prob = p * q * (p * P[2] +
      2 * (p * p) * P[3] +
      (1 / 2) * P[4] +
      p * P[5]);
  }
  if (left == 3 && right == 2) {
    prob = p * q * (q * P[2] +
      2 * (q * q) * P[3] +
      (1 / 2) * P[4] +
      q * P[5]);
  }
  if (left == 3 && right == 3) {
    prob = p * q * (2 * P[1] +
      P[2] +
      4 * p * q * P[3]);
  }
  return prob;
}


double calc_ll(double di,
               double l,
               double r,
               double t,
               int pop_size,
               double freq_ancestor_1,
               bool condition) {
  if (di < 0)
    return(-1e20);

  std::vector< double > seven_states = single_state(t, pop_size, di);

  std::vector< double > probs(3);
  double sum_prob = 0.0;
  for(int i = 0; i < 3; ++i) {
    probs[i] = get_prob_from_matrix(l,
                                    i,
                                    freq_ancestor_1,
                                    seven_states);
    sum_prob += probs[i];
  }

  double focal_prob = probs[r];
  if (condition == true)
    focal_prob *= 1.0 / sum_prob;

  return(log(focal_prob));
}

double chromosome::calculate_likelihood(double t, int pop_size, double freq_ancestor_1) {

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

 double ll = calc_ll(di, l, r, t, pop_size, freq_ancestor_1, false);


  for(int i = 2; i < distances.size(); ++i) {
    di = distances[i];
    l = states[i - 1];
    r = states[i];
    ll += calc_ll(di, l, r, t, pop_size, freq_ancestor_1, true);
  }
  return(ll);
}

