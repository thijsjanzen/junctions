#' @keywords internal
calc_ll_haploid_di <- function(info_vector, N, freq_ancestor_1, t) { # nolint
  di <- info_vector[[1]]
  left <- info_vector[[2]]
  right <- info_vector[[3]]

  H_0 <- 2 * freq_ancestor_1 * (1 - freq_ancestor_1) # nolint

  a1 <- H_0 * 2 * N / (2 * N + 1 / di) # nolint
  a2 <- 1 - (1 - di - 1 / (2 * N)) ^ t
  prob <- a1 * a2

  if (left == right) {
    prob <- 1 - prob
  }

  return(log(prob))
}

#' @keywords internal
calc_ll_haploid <- function(chrom_matrix,
                             N, # nolint
                             freq_ancestor_1,
                             t) {
  di <- c(diff(chrom_matrix[, 2]))
  to_analyze <- cbind(di,
                      chrom_matrix[1:(length(chrom_matrix[, 1]) - 1), 3],
                      chrom_matrix[2:length(chrom_matrix[, 1]), 3])

  ll <- apply(to_analyze, 1, calc_ll_haploid_di, N, freq_ancestor_1, t)
  return(sum(ll))
}

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
#' @return The number of generations passed since the onset of hybridization
#' @export
estimate_time_haploid <- function(ancestry_matrix,
                                  N = 1000, # nolint
                                  freq_ancestor_1 = 0.5,
                                  lower_lim = 2,
                                  upper_lim = 1000,
                                  verbose = FALSE) {

  if (!is.matrix(ancestry_matrix)) {
    ancestry_matrix <- as.matrix(ancestry_matrix)
  }

  chrom <- unique(ancestry_matrix[, 1])
  num_chrom <- length(chrom)

  calc_joint_ll <- function(t) {
    ll <- rep(0, num_chrom)
    for (i in seq_along(chrom)) {
      focal_chrom <- subset(ancestry_matrix, ancestry_matrix[, 1] == chrom[i])
      ll[i] <- calc_ll_haploid(focal_chrom, N, freq_ancestor_1, t)
    }
    if (verbose) cat(t, -sum(ll), "\n")
    return(-sum(ll))
  }

  fitted <- stats::optimize(calc_joint_ll, interval = c(lower_lim, upper_lim))
  if (fitted$minimum >= 0.9 * upper_lim) {
    warning("estimated time is close to the upper limit of time inference\n",
            "consider adjusting the upper limit to improve accuracy\n")
  }

  return(list(time = fitted$minimum,
              loglikelihood = -fitted$objective))
}
