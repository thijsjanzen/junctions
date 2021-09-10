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

#' log likelihood of the time since admixture for a haploid genome
#' @description log likelihood of the time since admixture for a set of single
#' chromosomes (for ex. in Yeast).
#' @param ancestry_matrix matrix with 3 columns, column 1 = chromosome,
#' column 2 = location in Morgan, column 3 = ancestry.
#' @param N Population Size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param t time since admixture
#' @return loglikelihood
#' @export
log_likelihood_haploid <- function(ancestry_matrix,
                                  N = 1000, # nolint
                                  freq_ancestor_1 = 0.5,
                                  t = 2) {
  ancestry_matrix <- as.matrix(ancestry_matrix)
  chrom <- unique(ancestry_matrix[, 1])
  num_chrom <- length(chrom)

  if (length(t) == 1) {
    ll <- rep(0, num_chrom)
    for (i in seq_along(chrom)) {
      focal_chrom <- subset(ancestry_matrix, ancestry_matrix[, 1] == chrom[i])
      ll[i] <- calc_ll_haploid(focal_chrom, N, freq_ancestor_1, t)
    }
    return(sum(ll))
  }
  if (length(t) > 1) {
    output_ll <- rep(NA, length(t))
    for (i in seq_along(t)) {
      ll <- rep(0, num_chrom)
      for (c in seq_along(chrom)) {
        focal_chrom <- subset(ancestry_matrix, ancestry_matrix[, 1] == chrom[c])
        ll[c] <- calc_ll_haploid(focal_chrom, N, freq_ancestor_1, t[i])
      }
      output_ll[i] <- sum(ll)
    }
    return(output_ll)
  }
}
