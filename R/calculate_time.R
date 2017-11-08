estimate_time <- function(J = NA, 
                          N = Inf, 
                          R = Inf, 
                          H_0 = 0.5, 
                          C = 1) {
  if (is.na(J)) {
    stop("ERROR! did you forget to provide J?\n")
  }

  if (is.infinite(N) && is.infinite(R)) {
    # following equation 1
    t <- J / (H_0 * C)
    return(t)
  }

  K <- junctions::calc_k(N, R, H_0, C)

  u <- 1 - 1 / (2 * N) - C / R

  t <- log(1 - J / K) / (log(u))
  return(t)
}