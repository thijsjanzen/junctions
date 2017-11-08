calculate_mat <- function(N = Inf, 
                          R = Inf, 
                          H_0 = 0.5, 
                          C = 1) {
  if (is.infinite(N) && is.infinite(R)) {
    stop("can not estimate MAT for both N and R infinite")
  }

  K <- junctions::calc_k(N, R, H_0, C)

  u <- 1 - 1 / (2 * N) - C / R
  mat <- log(1 / K) / log(u)
  return(mat)
}