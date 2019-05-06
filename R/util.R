#' @keywords internal
get_states <- function(local_anc_matrix) {
  # local anc is a matrix with 2 columns
  # column 1: ancestry chrom 1
  # column 2: ancestry chrom 2
  # we now generate a new vector, with the following coding:
  # 0 = (0, 0), homozygous for ancestor 0
  # 1 = (1, 1), homozygous for ancestor 1
  # 2 = (1, 2) or (2, 1), heterozygous

  local_anc <- rep(2, length(local_anc_matrix[,1]))
  homo_0 <- which(local_anc_matrix[,1] == 0 &
                    (local_anc_matrix[,1] == local_anc_matrix[,2]))
  homo_1 <- which(local_anc_matrix[,1] == 1 &
                    (local_anc_matrix[,1] == local_anc_matrix[,2]))
  local_anc[homo_0] <- 0
  local_anc[homo_1] <- 1

  all_states <- rep(NA, length(local_anc) - 1)
  for(i in 2:length(local_anc)) {
    left <- local_anc[i-1]
    right <- local_anc[i]

    if(left == 0 && right == 1) all_states[i-1] <- 1
    if(left == 1 && right == 0) all_states[i-1] <- 2
    if(left == 1 && right == 1) all_states[i-1] <- 3
    if(left == 0 && right == 0) all_states[i-1] <- 4
    if(left == 2 && right == 1) all_states[i-1] <- 5
    if(left == 2 && right == 0) all_states[i-1] <- 6
    if(left == 0 && right == 2) all_states[i-1] <- 7
    if(left == 1 && right == 2) all_states[i-1] <- 8
    if(left == 2 && right == 2) all_states[i-1] <- 9
  }
  return(all_states)
}

#' @keywords internal
single_state <- function(t, N, d) {
  trans_matrix = matrix(0,7,7)
  trans_matrix[1, ] = c(1 -1/(2*N) - 2*d , 2*d, 0, 0, 0, 1/(2*N), 0)
  trans_matrix[2, ] = c(1/(2*N), 1 - 3*1/(2*N) - d, d, 2*1/(2*N), 0, 0, 0)
  trans_matrix[3, ] = c(0, 2*1/(2*N), 1 - 4*1/(2*N), 0, 2 *1/(2*N), 0, 0)
  trans_matrix[4, ] = c(0, 0, 0, 1 - 1/(2*N) - d, d, 1/(2*N), 0)
  trans_matrix[5, ] = c(0, 0, 0, 2 * 1/(2*N), 1 - 3*1/(2*N), 0, 1/(2*N))
  trans_matrix[6, ] = c(0, 0, 0, 0, 0, 1 - d, d)
  trans_matrix[7, ] = c(0 ,0, 0, 0, 0, 1/(2*N), 1 - 1/(2*N))

  initial_state = c(1, 0, 0, 0, 0, 0, 0)

  output_state <- initial_state %*% expm::`%^%`(trans_matrix, t)
  return(output_state)
}

#' @keywords internal
get_expectation_O_state <- function(P, p, focal_state) {
  q <- 1-p
  cond_prob <- 1

  if(focal_state == 1) cond_prob = p*q*(p*q*P[3] + q*P[5] + P[7])

  if(focal_state == 2) cond_prob = p*q*(p*q*P[,3] + q*P[5] + P[7])

  if(focal_state == 3) cond_prob = (q^2)*( P[1] + P[4] + P[7]) +
                                   (q^3)*(P[2] + P[5]) +
                                   (q^4)*P[,3] + q*P[,6]


  if(focal_state == 4) cond_prob = (p^2)*( P[1] + P[,4] + P[,7]) +
                                   (p^3)*(P[2] + P[5]) +
                                   (p^4)*P[,3] + p*P[,6]

  if(focal_state == 5) cond_prob = p*q*(p*P[2] +
                                   2*(p^2)*P[3] +
                                   (1/2)*P[4] + p*P[5])

  if(focal_state == 6) cond_prob = p*q*(q*P[2] + 2*(q^2)*P[3] +
                                   (1/2)*P[4] + q*P[5])

  if(focal_state == 7) cond_prob = p*q*(q*P[2] + 2*(q^2)*P[3] +
                                   (1/2)*P[4] + q*P[5])

  if(focal_state == 8) cond_prob = p*q*(p*P[2] + 2*(p^2)*P[3] +
                                   (1/2)*P[4] + p*P[5])

  if(focal_state == 9) cond_prob = p*q*(2*P[1] + P[2] + 2*p*q*P[3])

  return(log(cond_prob))
}