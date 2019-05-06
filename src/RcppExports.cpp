// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// sim_fin_chrom
List sim_fin_chrom(int pop_size, double freq_ancestor_1, int run_time, double size_in_Morgan, int seed, int R);
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
List sim_inf_chrom(int pop_size, double freq_ancestor_1, int run_time, double size_in_Morgan, int markers, int seed);
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
List sim_phased_unphased_cpp(int pop_size, double freq_ancestor_1, int total_runtime, double size_in_morgan, int number_of_markers, NumericVector time_points, int seed);
RcppExport SEXP _junctions_sim_phased_unphased_cpp(SEXP pop_sizeSEXP, SEXP freq_ancestor_1SEXP, SEXP total_runtimeSEXP, SEXP size_in_morganSEXP, SEXP number_of_markersSEXP, SEXP time_pointsSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type freq_ancestor_1(freq_ancestor_1SEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type size_in_morgan(size_in_morganSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_markers(number_of_markersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type time_points(time_pointsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_phased_unphased_cpp(pop_size, freq_ancestor_1, total_runtime, size_in_morgan, number_of_markers, time_points, seed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_junctions_sim_fin_chrom", (DL_FUNC) &_junctions_sim_fin_chrom, 6},
    {"_junctions_sim_inf_chrom", (DL_FUNC) &_junctions_sim_inf_chrom, 6},
    {"_junctions_sim_phased_unphased_cpp", (DL_FUNC) &_junctions_sim_phased_unphased_cpp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_junctions(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
