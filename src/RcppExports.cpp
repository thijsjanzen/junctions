// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// simulate_backcrossing_cpp
Rcpp::List simulate_backcrossing_cpp(int pop_size, double freq_ancestor_1, int total_runtime, double size_in_morgan, int number_of_markers, Rcpp::NumericVector time_points, int seed);
RcppExport SEXP _junctions_simulate_backcrossing_cpp(SEXP pop_sizeSEXP, SEXP freq_ancestor_1SEXP, SEXP total_runtimeSEXP, SEXP size_in_morganSEXP, SEXP number_of_markersSEXP, SEXP time_pointsSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type freq_ancestor_1(freq_ancestor_1SEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type size_in_morgan(size_in_morganSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_markers(number_of_markersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type time_points(time_pointsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_backcrossing_cpp(pop_size, freq_ancestor_1, total_runtime, size_in_morgan, number_of_markers, time_points, seed));
    return rcpp_result_gen;
END_RCPP
}
// estimate_time_cpp
Rcpp::List estimate_time_cpp(const Rcpp::NumericMatrix& local_anc_matrix, const Rcpp::NumericVector& locations, int pop_size, double freq_ancestor_1, int lower_lim, int upper_lim, bool verbose, bool phased, int num_threads);
RcppExport SEXP _junctions_estimate_time_cpp(SEXP local_anc_matrixSEXP, SEXP locationsSEXP, SEXP pop_sizeSEXP, SEXP freq_ancestor_1SEXP, SEXP lower_limSEXP, SEXP upper_limSEXP, SEXP verboseSEXP, SEXP phasedSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type local_anc_matrix(local_anc_matrixSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type locations(locationsSEXP);
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type freq_ancestor_1(freq_ancestor_1SEXP);
    Rcpp::traits::input_parameter< int >::type lower_lim(lower_limSEXP);
    Rcpp::traits::input_parameter< int >::type upper_lim(upper_limSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type phased(phasedSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_time_cpp(local_anc_matrix, locations, pop_size, freq_ancestor_1, lower_lim, upper_lim, verbose, phased, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// loglikelihood_unphased_cpp
double loglikelihood_unphased_cpp(const Rcpp::NumericMatrix& local_anc_matrix, const Rcpp::NumericVector& locations, int pop_size, double freq_ancestor_1, double t, bool phased, bool verbose, int num_threads);
RcppExport SEXP _junctions_loglikelihood_unphased_cpp(SEXP local_anc_matrixSEXP, SEXP locationsSEXP, SEXP pop_sizeSEXP, SEXP freq_ancestor_1SEXP, SEXP tSEXP, SEXP phasedSEXP, SEXP verboseSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type local_anc_matrix(local_anc_matrixSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type locations(locationsSEXP);
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type freq_ancestor_1(freq_ancestor_1SEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< bool >::type phased(phasedSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikelihood_unphased_cpp(local_anc_matrix, locations, pop_size, freq_ancestor_1, t, phased, verbose, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// sim_fin_chrom
Rcpp::List sim_fin_chrom(int pop_size, double freq_ancestor_1, int run_time, double size_in_Morgan, int seed, int R);
RcppExport SEXP _junctions_sim_fin_chrom(SEXP pop_sizeSEXP, SEXP freq_ancestor_1SEXP, SEXP run_timeSEXP, SEXP size_in_MorganSEXP, SEXP seedSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type freq_ancestor_1(freq_ancestor_1SEXP);
    Rcpp::traits::input_parameter< int >::type run_time(run_timeSEXP);
    Rcpp::traits::input_parameter< double >::type size_in_Morgan(size_in_MorganSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_fin_chrom(pop_size, freq_ancestor_1, run_time, size_in_Morgan, seed, R));
    return rcpp_result_gen;
END_RCPP
}
// sim_inf_chrom
Rcpp::List sim_inf_chrom(int pop_size, double freq_ancestor_1, int run_time, double size_in_Morgan, int markers, int seed);
RcppExport SEXP _junctions_sim_inf_chrom(SEXP pop_sizeSEXP, SEXP freq_ancestor_1SEXP, SEXP run_timeSEXP, SEXP size_in_MorganSEXP, SEXP markersSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type freq_ancestor_1(freq_ancestor_1SEXP);
    Rcpp::traits::input_parameter< int >::type run_time(run_timeSEXP);
    Rcpp::traits::input_parameter< double >::type size_in_Morgan(size_in_MorganSEXP);
    Rcpp::traits::input_parameter< int >::type markers(markersSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_inf_chrom(pop_size, freq_ancestor_1, run_time, size_in_Morgan, markers, seed));
    return rcpp_result_gen;
END_RCPP
}
// sim_phased_unphased_cpp
Rcpp::List sim_phased_unphased_cpp(int pop_size, double freq_ancestor_1, int total_runtime, double size_in_morgan, Rcpp::NumericVector markers, Rcpp::NumericVector time_points, bool verbose, bool record_true_junctions, int num_indiv_sampled);
RcppExport SEXP _junctions_sim_phased_unphased_cpp(SEXP pop_sizeSEXP, SEXP freq_ancestor_1SEXP, SEXP total_runtimeSEXP, SEXP size_in_morganSEXP, SEXP markersSEXP, SEXP time_pointsSEXP, SEXP verboseSEXP, SEXP record_true_junctionsSEXP, SEXP num_indiv_sampledSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type freq_ancestor_1(freq_ancestor_1SEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type size_in_morgan(size_in_morganSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type markers(markersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type time_points(time_pointsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type record_true_junctions(record_true_junctionsSEXP);
    Rcpp::traits::input_parameter< int >::type num_indiv_sampled(num_indiv_sampledSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_phased_unphased_cpp(pop_size, freq_ancestor_1, total_runtime, size_in_morgan, markers, time_points, verbose, record_true_junctions, num_indiv_sampled));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_junctions_simulate_backcrossing_cpp", (DL_FUNC) &_junctions_simulate_backcrossing_cpp, 7},
    {"_junctions_estimate_time_cpp", (DL_FUNC) &_junctions_estimate_time_cpp, 9},
    {"_junctions_loglikelihood_unphased_cpp", (DL_FUNC) &_junctions_loglikelihood_unphased_cpp, 8},
    {"_junctions_sim_fin_chrom", (DL_FUNC) &_junctions_sim_fin_chrom, 6},
    {"_junctions_sim_inf_chrom", (DL_FUNC) &_junctions_sim_inf_chrom, 6},
    {"_junctions_sim_phased_unphased_cpp", (DL_FUNC) &_junctions_sim_phased_unphased_cpp, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_junctions(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
