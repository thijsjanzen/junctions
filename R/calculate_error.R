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

  error_t <- function(x) {
    if (is.na(x)) {
      return(NA)
    }

    K <- junctions::calc_k(N = N,      # nolint
                           R = R,      # nolint
                           H_0 = H_0,  # nolint
                           C = C)      # nolint

    u <- 1 - 1 / (2 * N) - C / R

    error <- log(u ^ x - 1 / K) / (log(u) * x) - 1
    # the flag relative determines whether we want the error
    # relative to K, or in absolute generations (relative = FALSE)
    if (relative) return(error)
    if (!relative) return(x * error)
  }

  output <- c()
  for (i in seq_along(t)) {
    output[i] <- error_t(t[i])
  }
  return(output)
}
