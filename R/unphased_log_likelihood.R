#' calculate the log likelihood of observing the data
#' @description Calculates the log likelihood of observing the unphased data, given the population size, initial heterozygosity and time since admixture
#' @param local_anc_matrix Local_anc can be provided as either matrix with two columns, where the first column represents ancestry on chromosome 1, and the second column represents ancestry on chromosome 2. Ancestry labels used should be [0, 1], where 0 indicates the first ancestor, and 1 indicates the second ancestor. Alternatively, the user can provide a vector indicating whether at the specific marker, the focal individual is homozygous for the first ancestor (0), homozygous for the second ancestor (1) or heterozygous (2).
#' @param locations locations of the used markers (in Morgan)
#' @param pop_size population size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param t time since admixture
#' @return log likelihood
#' @export
unphased_log_likelihood <- function(local_anc_matrix,
                                    locations,
                                    pop_size,
                                    freq_ancestor_1 = 0.5,
                                    t) {

  distances <- diff(locations)

  local_anc <- local_anc_matrix
  if(is.matrix(local_anc_matrix)) {
    local_anc <- rep(3, length(local_anc_matrix[,1]))
    homo_0 <- which(local_anc_matrix[,1] == 0 &
                      (local_anc_matrix[,1] == local_anc_matrix[,2]))
    homo_1 <- which(local_anc_matrix[,1] == 1 &
                      (local_anc_matrix[,1] == local_anc_matrix[,2]))
    local_anc[homo_0] <- 1
    local_anc[homo_1] <- 2
  }

  to_analyze <- cbind(distances,
                      local_anc[1:(length(local_anc)-1)],
                      local_anc[2:length(local_anc)]      )

  calc_ll <- function(local_time) {
    local_probs <- apply(to_analyze, 1, get_cond_prob_vector,
                         freq_ancestor_1,
                         pop_size,
                         local_time,
                         condition = TRUE)

    local_probs[1] <- get_cond_prob_vector(to_analyze[1,],
                                           freq_ancestor_1,
                                           pop_size,
                                           local_time,
                                           condition = FALSE)
    return(sum(local_probs))
  }

  if(length(t) == 1) {
    focal_ll <- calc_ll(t)

    return(focal_ll)
  }

  if(length(t) > 1) {
    output <- c()
    for(i in seq_along(t)) {
      focal_ll[i] <- calc_ll(t)
    }
    return(output)
  }
}