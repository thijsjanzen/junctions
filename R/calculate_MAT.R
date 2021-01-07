#' Function to calculate the maximum accurate time
#' @description Function that calculates the maximum time after hybridization
#' after which the number of junctions can still be reliably used to estimate
#' the onset of hybridization. This is following equation 15 in
#' Janzen et al. 2018.
#' @param N Population Size
#' @param R Number of genetic markers
#' @param H_0 Frequency of heterozygosity at t = 0
#' @param C Mean number of crossovers per meiosis (e.g. size in Morgan of
#'          the chromosome)
#' @return The maximum accurate time
#' @examples
#' calculate_mat(N = Inf, R = 1000, H_0 = 0.5, C = 1)
#' @keywords analytic time error
#' @export
calculate_mat <- function(N = Inf,   # nolint
                          R = Inf,   # nolint
                          H_0 = 0.5, # nolint
                          C = 1) {   # nolint
  if (is.infinite(N) && is.infinite(R)) {
    stop("can not estimate MAT for both N and R infinite")
  }

  K <- junctions::calc_k(N, R, H_0, C)  # nolint

  u <- 1 - 1 / (2 * N) - C / R
  mat <- log(1 / K) / log(u)
  return(mat)
}
