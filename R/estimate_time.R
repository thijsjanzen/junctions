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
  if (is.na(J)) {
    stop("ERROR! did you forget to provide J?\n")
  }

  if (is.infinite(N) && is.infinite(R)) {
    # following equation 1
    t <- J / (H_0 * C)
    return(t)
  }

  K <- junctions::calc_k(N, R, H_0, C)    # nolint

  u <- 1 - 1 / (2 * N) - C / R

  t <- log(1 - J / K) / (log(u))
  return(t)
}
