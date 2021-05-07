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
