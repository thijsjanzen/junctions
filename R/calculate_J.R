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
    jt <- H_0 * C * t
    return(jt)
  }

  K <- junctions::calc_k(N, R, H_0, C)

  jt <- K - K * (1 - H_0 * C / K) ^ t
  return(jt)
}

number_of_junctions_di <- function(N = Inf,
                                   H_0 = 0.5,
                                   C = 1,
                                   t = 100,
                                   di = 1e-6) {
 local_k <- 2 * N * C * H_0 / (2 * N * C + C / di)
 if(is.infinite(N)) local_k <- H_0 * C

 local_j <- local_k - local_k * (1 - 1 / (2*N) - di)^t
  return(local_j)
}


number_of_junctions_markers <- function(N = Inf,
                                        H_0 = 0.5,
                                        C = 1,
                                        t = 100,
                                        marker_distribution = NA) {

  if(length(marker_distribution) < 2) {
    stop("not provided a vector containing sufficient marker locations")
  }

  number_of_junctions <- 0
  for(i in 2:length(marker_distribution)) {
    di <- marker_distribution[i] - marker_distribution[i-1]
    expected_junctions <- number_of_junctions_di(N, H_0, C, t, di)
    number_of_junctions <- number_of_junctions + expected_junctions
  }
  return(number_of_junctions)
}


