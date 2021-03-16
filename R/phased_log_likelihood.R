#' calculate the log likelihood of observing the data
#' @description Calculates the log likelihood of observing the phased data,
#' given the population size, initial heterozygosity and time since admixture
#' @param local_anc_matrix Local_anc can be provided as either matrix with two
#' columns, where the first column represents ancestry on chromosome 1, and the
#' second column represents ancestry on chromosome 2. Ancestry labels used
#' should be [0, 1], where 0 indicates the first ancestor, and 1 indicates the
#' second ancestor. Alternatively, the user can provide a vector indicating
#' whether at the specific marker, the focal individual is homozygous for the
#' first ancestor (0), homozygous for the second ancestor (1) or
#' heterozygous (2).
#' @param locations locations of the used markers (in Morgan)
#' @param pop_size population size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param t time since admixture
#' @param num_threads number of threads, default is one thread. Set to -1 to
#'  use all available threads.
#' @return log likelihood
#' @export
loglikelihood_phased <- function(local_anc_matrix,
                                 locations,
                                 pop_size,
                                 freq_ancestor_1 = 0.5,
                                 t,
                                 num_threads = 1) {

  distances <- diff(locations)

  local_anc <- local_anc_matrix
  if (is.matrix(local_anc_matrix)) {
    local_anc <- rep(-1, length(local_anc_matrix[, 1]))
    homo_1 <- which(local_anc_matrix[, 1] == 0 &
                      local_anc_matrix[, 2] == 0)
    homo_2 <- which(local_anc_matrix[, 1] == 1 &
                      local_anc_matrix[, 2] == 1)

    het_1 <- which(local_anc_matrix[, 1] == 0 &
                     local_anc_matrix[, 2] == 1)
    het_2 <- which(local_anc_matrix[, 1] == 1 &
                     local_anc_matrix[, 2] == 0)

    local_anc[homo_1] <- 1
    local_anc[homo_2] <- 2
    local_anc[het_1]  <- 3
    local_anc[het_2]  <- 4
  }

  labels <- (unique(local_anc))
  if (sum(labels %in% 1:4) != length(labels)) {
    stop("local ancestry labels should be [1, 2, 3, 4] for homozygous P,
         homozygous Q, heterozygous PQ and heterozygous QP \n")
  }

  calc_ll <- function(params) {
    loglikelihood_unphased_cpp(local_anc_matrix = local_anc,
                               locations = locations,
                               pop_size = pop_size,
                               freq_ancestor_1 = freq_ancestor_1,
                               t = params[[1]],
                               phased = TRUE,
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