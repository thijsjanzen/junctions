#' calculate the log likelihood of observing the data
#' @description Calculates the log likelihood of observing the unphased data, given the population size, initial heterozygosity and time since admixture
#' @param local_anc matrix with two columns, where the first column represents ancestry on chromosome 1, and the second column represents ancestry on chromosome 2. Ancestry labels used should be [0, 1], where 0 indicates the first ancestor, and 1 indicates the second ancestor.
#' @param locations locations of the used markers (in Morgan)
#' @param pop_size population size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param t time since admixture
#' @return log likelihood
#' @export
unphased_log_likelihood <- function(local_anc,
                                    locations,
                                    pop_size,
                                    freq_ancestor_1 = 0.5,
                                    t) {

  distances <- diff(locations)
  local_states <- get_states(local_anc)

  calc_ll_single_state <- function(state, di,
                                   local_time,
                                   local_pop_size,
                                   local_p) {
    seven_states <- single_state(t = local_time, N = local_pop_size, d = di)
    focal_prob <- get_expectation_O_state(P = seven_states,
                                          p = local_p,
                                          focal_state = state)
    return(focal_prob)
  }

  local_probs <- mapply(calc_ll_single_state,
                        local_states,
                        distances,
                        t,
                        pop_size,
                        freq_ancestor_1)

  return(sum(local_probs))
}