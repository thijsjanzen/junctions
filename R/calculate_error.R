time_error <- function(J = NA,
                       N = Inf,
                       R = Inf,
                       H_0 = 0.5,
                       C = 1,
                       t = 1,
                       relative = TRUE) {
  # the flag relative determines whether we want the error
  # relative to K, or in absolute generations (relative = FALSE)
  K <- junctions::calc_k(N, R, H_0, C)

  u <- 1 - 1 / (2 * N) - C / R

  error <- log(u ^ t - 1 / K) / (log(u) * t) - 1

  if ( relative) return(  error)
  if (!relative) return(t*error)
}
