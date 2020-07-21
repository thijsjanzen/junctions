#' @keywords internal
get_prob_from_matrix <- function(left, right, p, P) {   # nolint
  # unphased probabilities!
  # left = marker on the left hand side [1 = PP, 2 = QQ, 3 = PQ or QP]
  # right = marker on the right hand side [1 = PP, 2 = QQ, 3 = PQ or QP]
  # p = frequency ancestor 1
  # P = states vector P_t^i

  q <- 1 - p
  prob <- 0
  if (left == 1 && right == 1) {
    prob <- (p ^ 2) * (P[1] + P[4] + P[7]) +
      (p ^ 3) * (P[2] + P[5]) +
      (p ^ 4) * P[3] +
      p * P[6]
  }
  if (left == 1 && right == 2) {
    prob <- p * q * (p * q * P[3] +
                    (1 / 2) * P[5] +
                    P[7])
  }
  if (left == 1 && right == 3) {
      prob <- p * q * (p * P[2] +
                     2 * (p ^ 2) * P[3] +
                   (1 / 2) * P[4] +
                   p * P[5])
  }

  if (left == 2 && right == 1) {
    prob <- p * q * (p * q * P[3] +
                   (1 / 2) * P[5] +
                   P[7])
  }
  if (left == 2 && right == 2) {
    prob <- (q ^ 2) * (P[1] + P[4] + P[7]) +
      (q ^ 3) * (P[2] + P[5]) +
      (q ^ 4) * P[3] +
      q * P[6]
  }
  if (left == 2 && right == 3) {
    prob <- p * q * (q * P[2] +
                   2 * (q ^ 2) * P[3] +
                   (1 / 2) * P[4] +
                   q * P[5])
  }

  if (left == 3 && right == 1) {
    prob <- p * q * (p * P[2] +
                   2 * (p ^ 2) * P[3] +
                   (1 / 2) * P[4] +
                   p * P[5])
  }
  if (left == 3 && right == 2) {
    prob <- p * q * (q * P[2] +
                   2 * (q ^ 2) * P[3] +
                   (1 / 2) * P[4] +
                   q * P[5])
  }
  if (left == 3 && right == 3) {
    prob <- p * q * (2 * P[1] +
                   P[2] +
                   4 * p * q * P[3])
  }
  return(prob)
}

#' @keywords internal
get_cond_prob_vector <- function(info_vector,
                                 freq_ancestor_1,
                                 pop_size,
                                 local_time,
                                 condition = TRUE) {
  di <- info_vector[[1]]
  left <- info_vector[[2]]
  right <- info_vector[[3]]

  if (local_time < 1)
    return(-Inf)
  if (pop_size < 2)
    return(-Inf)
  if (di < 0)
    return(-Inf)
  if (freq_ancestor_1 <= 0)
    return(-Inf)
  if (freq_ancestor_1 >= 1)
    return(-Inf)

  seven_states <- single_state(local_time, N = pop_size, d = di)

  probs <- c()
  for (j in 1:3) {
    probs[j] <- get_prob_from_matrix(left = left,
                                     right = j,
                                     p = freq_ancestor_1,
                                     P = seven_states)
  }

  focal_prob <- probs[right]

  if (condition == TRUE) {
    rel_prob <- focal_prob / sum(probs)
    final_prob <- log(rel_prob)
    return(final_prob)
  } else {
    return(log(focal_prob))
  }
}


#' estimates the time since admixture, given unphased ancestry data.
#' @description Calculates the time since admixture, given unphased
#' ancestry data.
#' @param local_anc_matrix Local_anc can be provided as either matrix\
#' with two columns, where the first column represents ancestry on
#' chromosome 1, and the second column represents ancestry on chromosome
#' 2. Ancestry labels used should be [0, 1], where 0 indicates the first
#' ancestor, and 1 indicates the second ancestor. Alternatively, the user
#' can provide a vector containing lables indicating whether at the
#' specific marker, the focal individual is homozygous for the first
#' ancestor (label = 1), homozygous for the second ancestor (label = 2)
#' or heterozygous (label = 3).
#' @param locations locations of the used markers (in Morgan)
#' @param pop_size population size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param lower_lim lower limit of the optimization algorithm. Increase
#' if the expected admixture time is relatively ancient
#' @param upper_lim upper limit of hte optimization algorithm. If set too
#' large, recent admixture events can be overlooked - best to set as low
#' as possible.
#' @param optim_pop_size If TRUE, population size is also optimized.
#' Starting point of the optimizaton will then be on the given population
#' size, and half upper_lim.
#' @param verbose display intermediate output? Default = FALSE
#' @export
estimate_time_unphased <- function(local_anc_matrix,
                                   locations,
                                   pop_size,
                                   freq_ancestor_1,
                                   lower_lim = 2,
                                   upper_lim = 1000,
                                   optim_pop_size = FALSE,
                                   verbose = FALSE) {

  distances <- diff(locations)

  local_anc <- local_anc_matrix
  if (is.matrix(local_anc_matrix)) {
    local_anc <- rep(3, length(local_anc_matrix[, 1]))
    homo_1 <- which(local_anc_matrix[, 1] == 0 &
                    local_anc_matrix[, 2] == 0)
    homo_2 <- which(local_anc_matrix[, 1] == 1 &
                    local_anc_matrix[, 2] == 1)

    local_anc[homo_1] <- 1
    local_anc[homo_2] <- 2
  }


  labels <- unique(local_anc)
  if (sum(labels %in% 1:3) != length(labels)) {
    stop("local ancestry labels should be [1, 2, 3] for homozygous anc 1,
         homozygous anc 2 and heterozygous\n")
  }

  to_analyze <- cbind(distances,
                      local_anc[1:(length(local_anc) - 1)],
                      local_anc[2:length(local_anc)])


  if (optim_pop_size == FALSE) {

    calc_ll <- function(params) {
      if (params[[1]] < 1) {
        return(Inf)
      }

      local_probs <- c()
      local_probs <- apply(to_analyze, 1, get_cond_prob_vector,
                           freq_ancestor_1,
                           pop_size,
                           local_time = params[[1]],
                           condition = TRUE)

      local_probs[1] <- get_cond_prob_vector(to_analyze[1, ],
                                             freq_ancestor_1,
                                             pop_size,
                                             local_time = params[[1]],
                                             condition = FALSE)

      if (verbose) cat(params[[1]], -sum(local_probs), "\n")
      return(-sum(local_probs))
    }

    a1 <- stats::optimize(f = calc_ll, interval = c(lower_lim, upper_lim))
    return(a1)
  }
  if (optim_pop_size == TRUE) {
    calc_ll <- function(params) {

      local_probs <- apply(to_analyze, 1, get_cond_prob_vector,
                           freq_ancestor_1,
                           params[[2]],
                           local_time = params[[1]],
                           condition = TRUE)

      local_probs[1] <- get_cond_prob_vector(to_analyze[1, ],
                                             freq_ancestor_1,
                                             params[[2]],
                                             local_time = params[[1]],
                                             condition = FALSE)

      if (verbose) cat(params[[1]], params[[2]], -sum(local_probs), "\n")
      return(-sum(local_probs))
    }

    a1 <- subplex::subplex(par = c(lower_lim, upper_lim), fn = calc_ll)
    return(a1)
  }
}
