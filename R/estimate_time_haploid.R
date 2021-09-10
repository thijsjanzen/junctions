#' estimate time using likelihood for a single chromosome
#' @description Estimate the time since the onset of hybridization, for a
#' haploid genome
#' @param ancestry_matrix matrix with 3 columns, column 1 = chromosome,
#' column 2 = location in Morgan, column 3 = ancestry.
#' @param N Population Size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param lower_lim lower limit of the optimization algorithm. Increase if the
#' expected admixture time is relatively ancient
#' @param upper_lim upper limit of the optimization algorithm. If set too large,
#' recent admixture events can be overlooked - best to set as low as possible.
#' @param verbose return verbose output
#' @param use_cpp use cpp version
#' @return The number of generations passed since the onset of hybridization
#' @export
estimate_time_haploid <- function(ancestry_matrix,
                                  N = 1000, # nolint
                                  freq_ancestor_1 = 0.5,
                                  lower_lim = 2,
                                  upper_lim = 1000,
                                  verbose = FALSE,
                                  use_cpp = FALSE) {

  if (!is.matrix(ancestry_matrix)) {
    ancestry_matrix <- as.matrix(ancestry_matrix)
  }

  chrom <- unique(ancestry_matrix[, 1])
  num_chrom <- length(chrom)

  fitted <- estimate_time_haploid_cpp(ancestry_matrix,
                                          N,
                                          freq_ancestor_1,
                                          lower_lim,
                                          upper_lim,
                                          verbose,
                                          use_cpp)

  return(list(time = fitted$time,
                loglikelihood = fitted$likelihood))
}
