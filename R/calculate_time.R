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


estimate_time_markers <- function(J = NA,
                                  N = Inf,
                                  H_0 = 0.5,
                                  marker_distribution = NA) {

  normalized_marker_distribution <- marker_distribution

  to_fit <- function(params) {
    expected_j <-
      number_of_junctions_markers(N = N,
                                  H_0 = H_0,
                                  t = params[[1]],
                                  marker_distribution =
                                        normalized_marker_distribution
                                  )
    return(abs(expected_j - J))
  }

  upper_lim = 1e5
  if(J > 1e4) {
    upper_lim = J * 20
  }

  fitted <- optimize(to_fit, interval = c(2, upper_lim) )
  return(fitted$minimum)
}
