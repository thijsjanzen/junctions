# nolint start
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
single_state_inf <- function(t, d) {

  p1 <- (1-d)^t
  p2 <- 2*(1-d)^t - 2*(1-2*d)^t
  p3 <- 1 + (1-2*d)^t - (1-d)^t
  return(  c(p1,p2,p3, 0, 0, 0, 0) )
}



# nolint end