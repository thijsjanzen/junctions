#' Estimate the time since the onset of hybridization, using the observed
#' number of junctions, taking into account the distribution of markers on
#' a single chromosome
#' @description Estimate the time since the onset of hybridization, following
#' equation 1 in Janzen et al. unpublished
#' @param J The observed number of junctions
#' @param N Population Size
#' @param H_0 Frequency of heterozygosity at t = 0
#' @param marker_distribution A vector containing the position of all markers
#' in Morgan.
#' @param lower_lim lower limit of the optimization algorithm. Increase if the
#' expected admixture time is relatively ancient
#' @param upper_lim upper limit of the optimization algorithm. If set too large,
#' recent admixture events can be overlooked - best to set as low as possible.
#' @return The number of generations passed since the onset of hybridization
#' @examples
#' markers <- seq(from = 0, to = 1, length.out = 100)
#' J <- number_of_junctions_markers(N = 100, H_0 = 0.5, t = 200,
#' marker_distribution = markers)
#' estimate_time_one_chrom(J = J,
#'                         N = 100,
#'                         H_0 = 0.5,
#'                         marker_distribution = markers) #should be 200 again
#' @export
estimate_time_one_chrom <- function(J = NA,          # nolint
                                    N = Inf,         # nolint
                                    H_0 = 0.5,       # nolint
                                    marker_distribution = NA,
                                    lower_lim = 2,
                                    upper_lim = 1000) {

  if (length(marker_distribution) < 2) {
   stop("not enough markers provided")
  }

  to_fit <- function(params) {
    expected_j <-
      number_of_junctions_markers(N = N,
                                  H_0 = H_0,
                                  t = params[[1]],
                                  marker_distribution)
    return(abs(expected_j - J))
  }

  fitted <- stats::optimize(to_fit, interval = c(lower_lim, upper_lim))
  if (fitted$minimum >= 0.9 * upper_lim) {
    warning("estimated time is close to the upper limit of time inference\n",
     "consider adjusting the upper limit to improve accuracy\n")
  }
  return(fitted$minimum)
}
