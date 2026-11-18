#' Estimate the time since the onset of hybridization, using the number of
#' junctions
#' @description Estimate the time since the onset of hybridization, following
#' equation 14 in Janzen et al. 2018
#' @param J The observed number of junctions
#' @param N Population Size
#' @param R Number of genetic markers
#' @param H_0 Frequency of heterozygosity at t = 0
#' @param C Mean number of crossovers per meiosis (e.g. size in Morgan of the
#' chromosome)
#' @return The number of generations passed since the onset of hybridization
#' @examples
#' J <- number_of_junctions(N = 100, R = 1000, H_0 = 0.5, C = 1, t = 200)
#' estimate_time(J = J, N = 100, R = 1000, H_0 = 0.5, C = 1)
#' # should be 200 again
#' @export
estimate_time <- function(J = NA,          # nolint
                          N = Inf,         # nolint
                          R = Inf,         # nolint
                          H_0 = 0.5,       # nolint
                          C = 1) {         # nolint
  if (all(is.na(J))) {
    stop("did you forget to provide J?")
  }
  args <- recycle_args(J, N, R, H_0, C)
  J <- args[[1]]  # nolint
  N <- args[[2]]  # nolint
  R <- args[[3]]  # nolint
  H_0 <- args[[4]]  # nolint
  C <- args[[5]]  # nolint

  K <- junctions::calc_k(N, R, H_0, C)    # nolint
  both_inf <- is.infinite(N) & is.infinite(R)

  too_many <- !both_inf & !is.na(J) & J >= K
  if (any(too_many)) {
    warning("NA value generated because the number of junctions is ",
            "higher than the expected maximum for this chromosome")
  }

  u <- 1 - 1 / (2 * N) - C / R
  t <- suppressWarnings(log(1 - J / K) / log(u))
  t[too_many] <- NA

  # following equation 1
  t[both_inf] <- J[both_inf] / (H_0[both_inf] * C[both_inf])
  return(t)
}
