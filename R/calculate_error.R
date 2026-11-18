#' Estimate the error in the time estimate
#' @description Calculate the error in the estimate of the onset of
#' hybridization, following Equations 3 & 4 in the Supplementary information of
#' Janzen et al. 2018.
#' @param t Inferred time
#' @param N Population Size
#' @param R Number of genetic markers
#' @param H_0 Frequency of heterozygosity at t = 0
#' @param C Mean number of crossovers per meiosis (e.g. size in Morgan of
#' the chromosome)
#' @param relative Boolean flag, if TRUE: return the relative error, if FALSE:
#' return error in generations
#' @return Expected error in the time estimate
#' @export
time_error <- function(t = NA,
                       N = Inf,    # nolint
                       R = Inf,    # nolint
                       H_0 = 0.5,  # nolint
                       C = 1,      # nolint
                       relative = TRUE) {
  args <- recycle_args(t, N, R, H_0, C)
  t <- args[[1]]
  N <- args[[2]]  # nolint
  R <- args[[3]]  # nolint
  H_0 <- args[[4]]  # nolint
  C <- args[[5]]  # nolint

  K <- junctions::calc_k(N = N, R = R, H_0 = H_0, C = C)   # nolint
  u <- 1 - 1 / (2 * N) - C / R

  # the error is only defined for t below the maximum accurate time (MAT),
  # i.e. while u^t > 1 / K
  beyond_mat <- !is.na(t) & u ^ t <= 1 / K
  if (any(beyond_mat)) {
    warning("NA returned for time points beyond the maximum accurate time ",
            "(see calculate_mat), where the error can not be calculated")
  }

  error <- suppressWarnings(log(u ^ t - 1 / K) / (log(u) * t) - 1)
  error[beyond_mat] <- NA

  # the flag relative determines whether we want the error
  # relative to K, or in absolute generations (relative = FALSE)
  if (relative) return(error)
  return(t * error)
}
