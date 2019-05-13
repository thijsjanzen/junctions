#' estimates the time since admixture, given unphased ancestry data.
#' @description Calculates the time since admixture, given unphased ancestry data.
#' @param local_anc Local_anc can be provided as either matrix with two columns, where the first column represents ancestry on chromosome 1, and the second column represents ancestry on chromosome 2. Ancestry labels used should be [0, 1], where 0 indicates the first ancestor, and 1 indicates the second ancestor. Alternatively, the user can provide a vector indicating whether at the specific marker, the focal individual is homozygous for the first ancestor (0), homozygous for the second ancestor (1) or heterozygous (2).
#' @param locations locations of the used markers (in Morgan)
#' @param pop_size population size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param max_t maximum time to be considered for maximum likelihood optimization. Too large values might cause the maximum likelihood algorithm to deviate from the global maximum.
#' @param optim_pop_size If TRUE, population size is also optimized. Starting point of the optimizaton will then be on the given population size, and half the maximum time.
#' @param verbose display intermediate output? Default = FALSE
#' @export
estimate_time_phased <- function(local_anc,
                                 locations,
                                 pop_size,
                                 freq_ancestor_1,
                                 max_t,
                                 optim_pop_size = FALSE,
                                 verbose = FALSE) {

  distances <- diff(locations)

  local_states <- get_states_phased(local_anc)

  calc_ll_single_state <- function(state, di,
                                   local_time,
                                   pop_size,
                                   freq_ancestor_1) {
    seven_states <- single_state(local_time, N = pop_size, d = di)
    focal_prob <- get_expectation_O_state_phased(seven_states,
                                                 p = freq_ancestor_1,
                                                 state)
    return(focal_prob)
  }

  if(optim_pop_size == FALSE) {

    to_optim_no_popsize <- function(t) {
      local_probs <- mapply(calc_ll_single_state,
                            local_states, distances,
                            local_time = t,
                            pop_size, freq_ancestor_1)
      if(verbose) cat(t, -sum(local_probs), "\n")
      return(-sum(local_probs))
    }
    a1 <- stats::optimize(f = to_optim_no_popsize, interval = c(2, max_t))
    return(a1)
  } else {
    to_optim_popsize <- function(params) {
      if(params[[1]] < 2) return(-Inf)
      if(params[[2]] < 2) return(-Inf)

      local_probs <- mapply(calc_ll_single_state,
                            local_states, distances,
                            local_time = params[[1]],
                            params[[2]], freq_ancestor_1)
      if(verbose) cat("N",params[[1]],"t", params[[2]], -sum(local_probs), "\n")
      return(-sum(local_probs))
    }
    a1 <- stats::optim(par = c(pop_size, max_t / 2), fn = to_optim_popsize)
    return(a1)
  }
}