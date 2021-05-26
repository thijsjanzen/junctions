#include <unistd.h>
#include <vector>

#include <cmath>
#include <iostream>

#include <Rcpp.h>
#include <nloptrAPI.h>

#include <RcppParallel.h>

#include "estimate_time_unphased_cpp.h"
#include "Output.h"


namespace detail {
int num_threads = 1;
}

struct chromosome {
  std::vector< size_t > states;
  std::vector< double > distances;
  bool phased;
  bool verbose;

  chromosome(const Rcpp::NumericVector& anc_matrix,
             const Rcpp::NumericVector& loc,
             bool p,
             bool verb) : phased(p), verbose(verb) {

    states = std::vector<size_t>(anc_matrix.begin(), anc_matrix.end());
    for(int i = 0; i < loc.size(); ++i) {
      if (i > 0) distances.push_back(loc[i] - loc[i - 1]);
    }
  }

  chromosome(const std::vector< std::vector< int > >& anc_matrix,
             const std::vector<double>& loc,
             bool p,
             bool verb) : phased(p), verbose(verb) {

    if (anc_matrix.size() != loc.size()) {
      Rcpp::stop("anc_matrix.nrow != loc.size()");
    }

    states = std::vector<size_t>(loc.size(), 2);

    if (phased) {
      for(size_t i = 0; i < loc.size(); ++i) {
        if (i > 0) distances.push_back(loc[i] - loc[i - 1]);
        if (anc_matrix[i][0] == anc_matrix[i][1]) {
          states[i] = anc_matrix[i][0];     // [0, 0] (state 0) or [1, 1] (state 1)
        } else {
          states[i] = 2 + anc_matrix[i][0] ; // [9, 1] (state 2) or [1, 0] (state 3)
        }
      }
    } else {
      for(size_t i = 0; i < anc_matrix.size(); ++i) {

        if (i > 0) {
          distances.push_back(loc[i] - loc[i - 1]);

          if (loc[i] - loc[i - 1] < 0) {
             Rcpp::stop("no negative distances allowed");
          }
        }

        if (anc_matrix[i][0] == anc_matrix[i][1]) {
          states[i] = anc_matrix[i][0];
        }
      }
    }
             }
  double calculate_likelihood(double t,
                              int pop_size,
                              double freq_ancestor_1) const;
};

struct nlopt_f_data {

  nlopt_f_data(const std::vector< chromosome >& c,
               int pop_size,
               double freq_anc) : chromosomes(c), N(pop_size), p(freq_anc) {
  }

  const std::vector< chromosome > chromosomes;
  const int N;
  const double p;
};

double objective(unsigned int n, const double* x, double*, void* func_data)
{
  auto psd = reinterpret_cast<nlopt_f_data*>(func_data);
  std::vector<double> ll(psd->chromosomes.size());
  for(size_t i = 0; i < psd->chromosomes.size(); ++i) {
    ll[i] = psd->chromosomes[i].calculate_likelihood(x[0], psd->N, psd->p);
  }
  auto sum_ll = -std::accumulate(ll.begin(), ll.end(), 0.0);
  if (psd->chromosomes.begin()->verbose) {
    Rcpp::Rcout << x[0] << " " << -sum_ll << "\n";
  }
  return sum_ll;
}

std::vector< chromosome > create_chromosomes(const Rcpp::NumericMatrix& local_anc_matrix,
                                             const Rcpp::NumericVector& locations,
                                             bool phased,
                                             bool verbose) {

  std::vector< chromosome > output;
  int current_chrom = local_anc_matrix(0, 0);

  std::vector< std::vector< int > > chrom_matrix;

 std::vector<double> positions;

  for(int i = 0; i < local_anc_matrix.nrow(); ++i) {
    bool add_chrom = false;
    if (local_anc_matrix(i, 0) != current_chrom) add_chrom = true;
    if (i > 0) {
      if (locations[i] < locations[i - 1]) add_chrom = true;
    }

    if (add_chrom) {
      current_chrom = local_anc_matrix(i, 0);
      chromosome new_chrom(chrom_matrix, positions, phased, verbose);
      output.push_back(new_chrom);
      chrom_matrix.clear();
      positions.clear();
    }
    std::vector<int> add = {static_cast<int>(local_anc_matrix(i, 1)),
                            static_cast<int>(local_anc_matrix(i, 2))};
    chrom_matrix.push_back(add);
    positions.push_back(locations[i]);
  }
  chromosome new_chrom(chrom_matrix, positions, phased, verbose);
  output.push_back(new_chrom);

  return output;
}



// [[Rcpp::export]]
Rcpp::List estimate_time_cpp(const Rcpp::NumericMatrix& local_anc_matrix,
                             const Rcpp::NumericVector& locations,
                             int pop_size,
                             double freq_ancestor_1,
                             int lower_lim,
                             int upper_lim,
                             bool verbose,
                             bool phased,
                             int num_threads = 1) {
try {

  if (verbose) {
    Rcpp::Rcout << "welcome to estimate_time_cpp\n";
  }

  detail::num_threads = num_threads;

  if (local_anc_matrix.ncol() != 3) {
    Rcpp::stop("local ancestry matrix has to have 3 columns");
  }


  if (verbose) {
    Rcpp::Rcout << "starting create chromosomes\n";
  }
  std::vector< chromosome > chromosomes = create_chromosomes(local_anc_matrix,
                                                             locations,
                                                             phased,
                                                             verbose);

  if (verbose) {
    Rcpp::Rcout << "chromosomes read from data\n";
  }

  nlopt_f_data optim_data(chromosomes, pop_size, freq_ancestor_1);

  nlopt_opt opt = nlopt_create(NLOPT_LN_SBPLX, static_cast<unsigned int>(1));

  double llim[1] = {static_cast<double>(lower_lim)};
  double ulim[1] = {static_cast<double>(upper_lim)};

  nlopt_set_lower_bounds(opt, llim);
  nlopt_set_upper_bounds(opt, ulim);

  nlopt_set_min_objective(opt, objective, &optim_data);

  nlopt_set_xtol_rel(opt, 1e-1);
  std::vector<double> x = {10};
  double minf;
  if (verbose) {
    Rcpp::Rcout << "starting optimisation\n";
  }

  auto nloptresult = nlopt_optimize(opt, &(x[0]), &minf);

  if (nloptresult < 0) {
    Rcpp::Rcout << "failure to optimize!\n";
  }

  if (verbose) {
    Rcpp::Rcout << "done optimisation\n";
  }

  nlopt_destroy(opt);

  std::vector<double> output = {x[0], minf};

  return Rcpp::List::create(Rcpp::Named("time") = x[0],
                            Rcpp::Named("likelihood") = -1 * minf);
} catch(std::exception &ex) {
  forward_exception_to_r(ex);
} catch(...) {
  ::Rf_error("c++ exception (unknown reason)");
}
return NA_REAL;
}

// [[Rcpp::export]]
double loglikelihood_unphased_cpp(const Rcpp::NumericMatrix& local_anc_matrix,
                                  const Rcpp::NumericVector& locations,
                                  int pop_size,
                                  double freq_ancestor_1,
                                  double t,
                                  bool phased,
                                  bool verbose = false,
                                  int num_threads = 1) {

  detail::num_threads = num_threads;
  if (local_anc_matrix.ncol() != 3) {
    Rcpp::stop("local ancestry matrix has to have 3 columns");
  }

  std::vector< chromosome > chromosomes = create_chromosomes(local_anc_matrix,
                                                             locations,
                                                             phased,
                                                             verbose);

  std::vector< double > ll(chromosomes.size());
  for (size_t i = 0; i < chromosomes.size(); ++i) {
    ll[i] = chromosomes[i].calculate_likelihood(t, pop_size, freq_ancestor_1);
  }

  return std::accumulate(ll.begin(), ll.end(), 0.0);
}

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


  SqMx<7, double> m(trans_matrix);
  SqMx<7, double> m2 = m ^ t;
  std::vector< double > output = m2();

  return(output);
}


double get_prob_from_matrix_unphased_cpp(int left,
                                         int right,
                                         double p,
                                         const std::vector<double>& P) {
  // unphased probabilities!
  // left = marker on the left hand side [1 = PP, 2 = QQ, 3 = PQ or QP]
  // right = marker on the right hand side [1 = PP, 2 = QQ, 3 = PQ or QP]
  // p = frequency ancestor 1
  // P = states vector P_t^i

  left++; // to conform to R notation. code below is a close copy of R code
  right++;

  double q  = 1 - p;
  double prob = 0;

  if (left == 1 && right == 1) {
    prob = (p * p) * (P[0] + P[3] + P[6]) +
      pow(p, 3) * (P[1] + P[4]) +
      pow(p, 4) * P[2] +
      p * P[5];
  }
  if (left == 1 && right == 2) {
    prob = p * q * (p * q * P[2] +
      (1.0 / 2) * P[4] +
      P[6]);
  }
  if (left == 1 && right == 3) {
    prob = p * q * (p * P[1] +
      2 * (p * p) * P[2] +
      (1.0 / 2) * P[3] +
      p * P[4]);
  }

  if (left == 2 && right == 1) {
    prob = p * q * (p * q * P[2] +
      (1.0 / 2) * P[4] +
      P[6]);
  }
  if (left == 2 && right == 2) {
    prob = (q * q) * (P[0] + P[3] + P[6]) +
      pow(q, 3) * (P[1] + P[4]) +
      pow(q, 4) * P[2] +
      q * P[5];
  }

  if (left == 2 && right == 3) {
    prob = p * q * (q * P[1] +
      2 * (q * q) * P[2] +
      (1.0 / 2) * P[3] +
      q * P[4]);
  }

  if (left == 3 && right == 1) {
    prob = p * q * (p * P[1] +
      2 * (p * p) * P[2] +
      (1.0 / 2) * P[3] +
      p * P[4]);
  }
  if (left == 3 && right == 2) {
    prob = p * q * (q * P[1] +
      2 * (q * q) * P[2] +
      (1.0 / 2) * P[3] +
      q * P[4]);
  }
  if (left == 3 && right == 3) {
    prob = p * q * (2 * P[0] +
      P[1] +
      4 * p * q * P[2]);
  }
  return prob;
}

double get_prob_from_matrix_phased_cpp(int left,
                                       int right,
                                       double p,
                                       const std::vector<double>& P) {
  // unphased probabilities!
  // left = marker on the left hand side [1 = PP, 2 = QQ, 3 = PQ or QP]
  // right = marker on the right hand side [1 = PP, 2 = QQ, 3 = PQ or QP]
  // p = frequency ancestor 1
  // P = states vector P_t^i

  left++; // to conform to R notation. code below is a close copy of R code
  right++;

  double q  = 1 - p;
  double prob = 0;

  if (left == 1 && right == 1) {
    double  p_sq = p * p;
    prob =  p_sq * (P[0] + P[3] + P[6]) +
      (p_sq * p) * (P[1] + P[4]) +
      (p_sq * p_sq) * P[2] +
      p * P[5];
  }
  if (left == 1 && right == 2) {
    prob =  p * q * (p * q * P[2] +
      (1.0 / 2.0) * P[4] +
      P[6]);
  }
  if (left == 1 && right == 3) {
    prob =  (1.0 / 2.0) * (p * q) * (p * P[1] +
      2 * (p * p) * P[2] +
      (1.0 / 2.0) * P[3] +
      p * P[4]);
  }

  if (left == 1 && right == 4) {
    prob =  (1.0 / 2.0) * (p * q) * (p * P[1] +
      2 * (p * p) * P[2] +
      (1.0 / 2.0) * P[3] +
      p * P[4]);
  }


  if (left == 2 && right == 1) {
    prob =  p * q * (p * q * P[2] +
      (1.0 / 2.0) * P[4] +
      P[6]);
  }
  if (left == 2 && right == 2) {

    double q_sq = q * q;
    prob =  (q_sq) * (P[0] + P[3] + P[6]) +
      (q_sq * q) * (P[1] + P[4]) +
      (q_sq * q_sq) * P[2] +
      q * P[5];
  }

  if (left == 2 && right == 3) {
    prob =  (1.0 / 2.0) * p * q * (q * P[1] +
      2 * (q * q) * P[2] +
      (1.0 / 2.0) * P[3] +
      q * P[4]);
  }

  if (left == 2 && right == 4) {
    prob =  (1.0 / 2.0) * p * q * (q * P[1] +
      2 * (q * q) * P[2] +
      (1.0 / 2.0) * P[3] +
      q * P[4]);
  }

  if (left == 3 && right == 1) {
    prob =  (1.0 / 2.0) * p * q * (p * P[1] +
      2 * (p * p) * P[2] +
      (1.0 / 2.0) * P[3] +
      p * P[4]);
  }
  if (left == 3 && right == 2) {
    prob =  (1.0 / 2.0) * p * q * (q * P[1] +
      2 * (q * q) * P[2] +
      (1.0 / 2.0) * P[3] +
      q * P[4]);
  }
  if (left == 3 && right == 3) {
    prob =  p * q * (P[0] +
      (1.0 / 2.0) * P[1] +
      p * q * P[2]);
  }
  if (left == 3 && right == 4) {
    prob =  (p * p) * (q * q) * P[2];
  }

  if (left == 4 && right == 1) {
    prob =  (1.0 / 2.0) * p * q * (p * P[1] +
      2 * (p * p) * P[2] +
      (1.0 / 2.0) * P[3] +
      p * P[4]);
  }
  if (left == 4 && right == 2) {
    prob =  (1.0 / 2.0) * p * q * (q * P[1] +
      2 * (q * q) * P[2] +
      (1.0 / 2.0) * P[3] +
      q * P[4]);
  }
  if (left == 4 && right == 3) {
    prob =  (p * p) * (q * q) * P[2];
  }
  if (left == 4 && right == 4) {
    prob =  (1.0 / 2.0) * p * q * (2 * P[0] +
      P[1] +
      2 * p * q * P[2]);
  }

  return prob;
}

double calc_ll(double di,
               double l,
               double r,
               double t,
               int pop_size,
               double freq_ancestor_1,
               bool condition,
               bool phased) {
  if (di < 0) {
    Rcpp::Rcout << "di < 0\n";
    return(-1e20);
  }

  std::vector<  double > seven_states = single_state_cpp(t, pop_size, di);

  int num_dim = 3 + phased;

  std::vector<  double > probs(num_dim);
  double sum_prob = 0.0;
  for(int i = 0; i < num_dim; ++i) {
    if (phased) {
      probs[i] = get_prob_from_matrix_phased_cpp(l,
                                                 i,
                                                 freq_ancestor_1,
                                                 seven_states);
    } else {
      probs[i] = get_prob_from_matrix_unphased_cpp(l,
                                                   i,
                                                   freq_ancestor_1,
                                                   seven_states);
    }
    sum_prob += probs[i];
  }

  double focal_prob = probs[r];

  if (condition == true)
    focal_prob *= 1.0 / sum_prob;

  return(log(focal_prob));
}

double chromosome::calculate_likelihood(double t,
                                        int pop_size,
                                        double freq_ancestor_1) const {

  if (t < 1) {
    Rcpp::Rcout << "t < 1\n";
    return(-1e20);
  }
  if (pop_size < 2) {
    Rcpp::Rcout << "pop_size < 2\n";
    return(-1e20);
  }
  if (freq_ancestor_1 >= 1) {
    Rcpp::Rcout << "p >= 1\n";
    return(-1e20);
  }
  if (freq_ancestor_1 <= 0) {
    Rcpp::Rcout << "p <= 0\n";
    return(-1e20);
  }

  std::vector< double> ll(distances.size());
  ll[0] = calc_ll(distances[0],
                  states[0],
                        states[1], t, pop_size, freq_ancestor_1, false,
                        phased);

  if (detail::num_threads == 1) {
    for (size_t i = 0; i < distances.size(); ++i) {

      double di = distances[i];
      double l = states[i];
      double r = states[i + 1];
      ll[i] = calc_ll(di, l, r, t, pop_size, freq_ancestor_1, true,
                      phased);
    }
  } else {

    tbb::task_scheduler_init _tbb((detail::num_threads > 0) ? detail::num_threads : tbb::task_scheduler_init::automatic);

    tbb::parallel_for(
      tbb::blocked_range<unsigned>(1, distances.size()),
      [&](const tbb::blocked_range<unsigned>& r) {
        for (unsigned i = r.begin(); i < r.end(); ++i) {

          double di = distances[i];
          double l = states[i];
          double r = states[i + 1];
          ll[i] = calc_ll(di, l, r, t, pop_size, freq_ancestor_1, true,
                          phased);
        }
      }
    );
  }
  double answer = std::accumulate(ll.begin(), ll.end(), 0.0);
  return(answer);
}
