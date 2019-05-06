#' estimates the time since admixture, given unphased ancestry data.
#' @description Calculates the time since admixture, given unphased ancestry data.
#' @param local_anc matrix with two columns, where the first column represents ancestry on chromosome 1, and the second column represents ancestry on chromosome 2. Ancestry labels used should be [0, 1], where 0 indicates the first ancestor, and 1 indicates the second ancestor.
#' @param locations locations of the used markers (in Morgan)
#' @param pop_size population size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param max_t maximum time to be considered for maximum likelihood optimization. Too large values might cause the maximum likelihood algorithm to deviate from the global maximum.
#' @export
estimate_time_unphased <- function(local_anc,
                                   locations,
                                   pop_size,
                                   freq_ancestor_1,
                                   max_t) {

  distances <- diff(locations)
  local_states <- get_states(local_anc)

  calc_ll_single_state <- function(state, di,
                                   local_time,
                                   pop_size,
                                   freq_ancestor_1) {
    seven_states <- single_state(local_time, N = pop_size, d = di)
    focal_prob <- get_expectation_O_state(seven_states, p = freq_ancestor_1, state)
    return(focal_prob)
  }

  to_optim <- function(t) {
    local_probs <- mapply(calc_ll_single_state,
                          local_states, distances,
                          local_time = t,
                          pop_size, freq_ancestor_1)
    #cat(t, -sum(local_probs), "\n")
    return(-sum(local_probs))
  }
  a1 <- stats::optimize(f = to_optim, interval = c(2, max_t))
  return(a1)
}