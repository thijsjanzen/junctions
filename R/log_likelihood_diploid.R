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

  x <- as.data.frame(local_anc_matrix)
  if (ncol(x) != 4) {
    stop("local_anc_matrix must have four columns: chromosome, location, ",
         "ancestry chromosome 1, ancestry chromosome 2")
  }
  # re-use the diploid input checks, treating all data as one individual
  input <- prepare_diploid_input(cbind(1, x))
  m <- input$data

  RcppParallel::setThreadOptions(num_threads)

  vapply(t, function(focal_t) {
    loglikelihood_unphased_cpp(local_anc_matrix = m[, c(6, 4, 5),
                                                    drop = FALSE],
                               locations = m[, 3],
                               pop_size = pop_size,
                               freq_ancestor_1 = freq_ancestor_1,
                               t = focal_t,
                               phased = phased)
  }, numeric(1))
}
