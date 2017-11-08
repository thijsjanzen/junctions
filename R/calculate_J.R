number_of_junctions <- function(N = Inf,
                                R = Inf,
                                H_0 = 0.5,
                                C = 1,
                                t = 100) {

  if (is.infinite(N) && is.infinite(R)) {
    # If both N and R are infinite, R gives
    # numerical problems using equation 12
    # to calculate K, so instead we use
    # equation 1
    jt <- H_0 * C * t;
    return(jt)
  }

  K <- junctions::calc_k(N, R, H_0, C)

  jt <- K - K * (1 - H_0 * C / K) ^ t
  return(jt)
}
