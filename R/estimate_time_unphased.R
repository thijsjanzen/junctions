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
estimate_time_unphased <- function(local_anc,
                                   locations,
                                   pop_size,
                                   freq_ancestor_1,
                                   max_t,
                                   optim_pop_size = FALSE,
                                   verbose = FALSE) {

  distances <- diff(locations)

  local_states <- get_states(local_anc)

  calc_ll_single_state <- function(state, di,
                                   local_time,
                                   pop_size,
                                   freq_ancestor_1) {
    seven_states <- single_state(local_time, N = pop_size, d = di)
    focal_prob <- get_expectation_O_state(seven_states,
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

#' @keywords internal
get_cond_prob <- function(left,
                          right,
                          P,
                          freq_ancestor_1) {

  state_matrix <- matrix(1:9, nrow = 3, ncol = 3, byrow = T)
  focal_state <- state_matrix[left, right]
  focal_prob <- get_expectation_O_state(P,
                                        p = freq_ancestor_1,
                                        focal_state)
  focal_prob <- exp(focal_prob)

  norm_constant <- 0
  for(j in 1:3) {
    norm_constant <- norm_constant + exp(get_expectation_O_state(P,
                                                             p = freq_ancestor_1,
                                                             state_matrix[left, j]))
  }
  final_prob <- focal_prob / norm_constant
  return(log(final_prob))
}

#' @keywords internal
get_cond_prob_vector <- function(info_vector,
                                 freq_ancestor_1,
                                 pop_size,
                                 local_time) {
  di <- info_vector[1]
  left = info_vector[2]
  right = info_vector[3]

  seven_states <- single_state(local_time, N = pop_size, d = di)

  state_matrix <- matrix(1:9, nrow = 3, ncol = 3, byrow = T)

  log_probs <- c()
  for(j in 1:3) {
    log_probs[j] <- get_expectation_O_state(seven_states,
                                            p = freq_ancestor_1,
                                            state_matrix[left, j])
  }

  log_probs <- exp(log_probs)

  focal_prob <- log_probs[right]

  final_prob <- log( focal_prob / sum(log_probs) )
  #cat(di, left, right, log_probs, exp(final_prob), "\n")

  return(final_prob)
}


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
estimate_time_unphased2 <- function(local_anc_matrix,
                                   locations,
                                   pop_size,
                                   freq_ancestor_1,
                                   max_t,
                                   optim_pop_size = FALSE,
                                   verbose = FALSE) {

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
                      local_anc[2:length(local_anc)])


  calc_ll <- function(params) {

    local_probs <- apply(to_analyze, 1, get_cond_prob_vector,
                         freq_ancestor_1,
                         pop_size,
                         local_time = params[[1]])

    if(verbose) cat(params[[1]], -sum(local_probs), "\n")
    return(-sum(local_probs))
  }

  if(1 == 2) {
    found_ll <- c()
    focal_t <- c(1:50, seq(100, 1000, by = 100), seq(1000, 10000, by = 500))
    for(i in seq_along(focal_t)) {
      found_ll[i] <- calc_ll(focal_t[i])
    }
    plot(found_ll~focal_t, log = "x")
  }
  a1 <- stats::optimize(f = calc_ll, interval = c(2, max_t))
  return(a1)
}