# nolint start
#' @keywords internal
get_num_markers <- function(markers) {
  if (length(markers) == 1) {
    num_markers <- abs(markers[1])

    if (markers[1] < 0) { # evenly spaced markers
      di = 1.0 / (num_markers);
      markers <- seq(di, 1 - di, length.out = num_markers)
      return(markers)
    } else {
      markers <- stats::runif(n = num_markers)
      markers <- unique(markers)
      while (length(markers) < num_markers) {
        add <- stats::runif(num_markers - length(markers))
        markers <- unique(c(markers, add))
      }
      markers <- sort(markers)
      return(markers)
    }
  } else {
    return(markers)
  }
}

#' @keywords internal
single_state <- function(t, N, d) {

  trans_matrix <- matrix(0, 7, 7)
  trans_matrix[1, ] <- c(1 - 1 / (2*N) - 2 * d , 2 * d, 0, 0, 0, 1 / (2*N), 0)
  trans_matrix[2, ] <- c(1 / (2*N), 1 - 3 * 1 / (2*N) - d, d, 2 * 1 / (2*N), 0, 0, 0)
  trans_matrix[3, ] <- c(0, 2 * 1 / (2*N), 1 - 4 * 1/(2*N), 0, 2 * 1/(2*N), 0, 0)
  trans_matrix[4, ] <- c(0, 0, 0, 1 - 1 / (2*N) - d, d, 1/(2*N), 0)
  trans_matrix[5, ] <- c(0, 0, 0, 2 * 1 / (2*N), 1 - 3 * 1/(2*N), 0, 1 / (2*N))
  trans_matrix[6, ] <- c(0, 0, 0, 0, 0, 1 - d, d)
  trans_matrix[7, ] <- c(0 ,0, 0, 0, 0, 1 / (2*N),  1 - 1 / (2*N))

  initial_state <- c(1, 0, 0, 0, 0, 0, 0)

  output_state <- initial_state %*% expm::`%^%`(trans_matrix, t)
  return(output_state)
}

#' @keywords internal
get_index_and_prob_unphased <- function(info_vector,
                               freq_ancestor_1,
                               pop_size,
                               local_time) {
  di <- info_vector[[1]]
  left <- info_vector[[2]]
  right <- info_vector[[3]]

  seven_states <- single_state(local_time, N = pop_size, d = di)

  probs <- c()
  for (j in 1:3) {
    probs[j] <- get_prob_from_matrix(left = left,
                                     right = j,
                                     p = freq_ancestor_1,
                                     P = seven_states)
  }

  focal_prob <- probs[right]

  rel_prob <- focal_prob / sum(probs)
  final_prob <- log(rel_prob)

  return(c(left, right, final_prob))
}


#' @keywords internal
get_index_and_prob_phased <- function(info_vector,
                                        freq_ancestor_1,
                                        pop_size,
                                        local_time) {
  di <- info_vector[[1]]
  left <- info_vector[[2]]
  right <- info_vector[[3]]

  seven_states <- single_state(local_time, N = pop_size, d = di)

  probs <- c()
  for (j in 1:4) {
    probs[j] <- get_prob_from_matrix_phased(left = left,
                                            right = j,
                                            p = freq_ancestor_1,
                                            P = seven_states)
  }

  focal_prob <- probs[right]

  rel_prob <- focal_prob / sum(probs)
  final_prob <- log(rel_prob)

  return(c(left, right, final_prob))
}


#' @keywords internal
matrix_of_states_unphased <- function(local_anc_matrix,
           locations,
           pop_size,
           freq_ancestor_1 = 0.5,
           t) {

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


    labels <- (unique(local_anc))
    if (sum(labels %in% 1:3) != length(labels)) {
      stop("local ancestry labels should be [1, 2, 3] for homozygous anc 1,
         homozygous anc 2 and heterozygous\n")
    }

    to_analyze <- cbind(distances,
                        local_anc[1:(length(local_anc) - 1)],
                        local_anc[2:length(local_anc)])

    # we need a new matrix with:
    # pos, pos in matrix, value
    results <- apply(to_analyze, 1, get_index_and_prob_unphased,
                     freq_ancestor_1,
                     pop_size,
                     t)
    results <- t(results)

    ll_matrix <- matrix(0, nrow = 3, ncol = 3)
    for(i in 1:length(results[,1])) {
      index_left <- results[i, 1]
      index_right <- results[i, 2]
      val   <- results[i, 3]
      ll_matrix[index_left, index_right] <- ll_matrix[index_left, index_right] + val
    }

    return(ll_matrix)
}

#' @keywords internal
matrix_of_states_phased <- function(local_anc_matrix,
                                      locations,
                                      pop_size,
                                      freq_ancestor_1 = 0.5,
                                      t) {

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

  to_analyze <- cbind(distances,
                      local_anc[1:(length(local_anc) - 1)],
                      local_anc[2:length(local_anc)])

  # we need a new matrix with:
  # pos, pos in matrix, value
  results <- apply(to_analyze, 1, get_index_and_prob_phased,
                   freq_ancestor_1,
                   pop_size,
                   t)

  results <- t(results)

  ll_matrix <- matrix(0, nrow = 4, ncol = 4)
  for(i in 1:length(results[,1])) {
    index_left <- results[i, 1]
    index_right <- results[i, 2]
    val   <- results[i, 3]
    ll_matrix[index_left, index_right] <- ll_matrix[index_left, index_right] + val
  }

  return(ll_matrix)
}

#' @keywords internal
count_of_states_phased <- function(local_anc_matrix) {

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

  count_matrix <- matrix(0, nrow = 4, ncol = 4)
  for(i in 1:(length(local_anc) - 1)) {
    index_left <- local_anc[i]
    index_right <- local_anc[i + 1]
    count_matrix[index_left, index_right] <- count_matrix[index_left, index_right] + 1
  }

  return(count_matrix)
}


#' @keywords internal
count_of_states_unphased <- function(local_anc_matrix) {

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

  count_matrix <- matrix(0, nrow = 3, ncol = 3)
  for(i in 1:(length(local_anc) - 1)) {
    index_left <- local_anc[i]
    index_right <- local_anc[i + 1]
    count_matrix[index_left, index_right] <- count_matrix[index_left, index_right] + 1
  }

  return(count_matrix)
}

# nolint end