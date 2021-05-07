#' calculate the log likelihood of observing diploid ancestry data.
#' @description Calculates the log likelihood of observing the phased data,
#' given the population size, initial heterozygosity and time since admixture
#' @param local_anc_matrix a matrix with four columns: column 1) chromosome
#' indicator, 2) location of marker in Morgan on respective chromosome 3)
#' ancestry at chromosome 4) ancestry at chromosome 2.
#' @param pop_size population size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param t time since admixture
#' @param phased is the data phased or not? default is false.
#' @param num_threads number of threads, default is one thread. Set to -1 to
#'  use all available threads.
#' @return log likelihood
#' @export
log_likelihood_diploid <- function(local_anc_matrix,
                                  pop_size,
                                  freq_ancestor_1 = 0.5,
                                  t,
                                  phased = FALSE,
                                  num_threads = 1) {

  local_anc_matrix <- as.matrix(local_anc_matrix)
  locations <- local_anc_matrix[, 2]
  local_anc <- local_anc_matrix[, c(1, 3, 4)]

  calc_ll <- function(params) {
    loglikelihood_unphased_cpp(local_anc_matrix = local_anc,
                               locations = locations,
                               pop_size = pop_size,
                               freq_ancestor_1 = freq_ancestor_1,
                               t = params[[1]],
                               phased = phased,
                               num_threads = num_threads)
  }

  if (length(t) == 1) {
    focal_ll <- calc_ll(t)

    return(focal_ll)
  }

  if (length(t) > 1) {
    output <- c()
    for (i in seq_along(t)) {
      output[i] <- calc_ll(t[i])
    }
    return(output)
  }
}
