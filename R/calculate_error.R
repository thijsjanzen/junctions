#' Estimate the error in the time estimate
#' @description Calculate the error in the estimate of the onset of hybridisation, following Equations 3 & 4 in the Supplementary information of Janzen et al. 2018.
#' @param J The number of junctions at time t
#' @param N Population Size
#' @param R Number of genetic markers
#' @param H_0 Frequency of heterozygosity at t = 0
#' @param C Mean number of crossovers per meiosis (e.g. size in Morgan of the chromosome)
#' @param t Inferred time
#' @param relative Boolean flag, if TRUE: return the relative error, if FALSE: return error in generations
#' @return Expected error in the time estimate
#' @examples
#' time_error(J = 100, N = Inf, R = 1000, H_0 = 0.5, C = 1)
#' @keywords analytic time error
#' @export
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
