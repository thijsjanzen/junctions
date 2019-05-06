#' Estimate the time since the onset of hybridization, using the observed number of junctions, taking into account the distribution of markers
#' @description Estimate the time since the onset of hybridization, following equation 1 in Janzen et al. unpublished
#' @param J The observed number of junctions
#' @param N Population Size
#' @param H_0 Frequency of heterozygosity at t = 0
#' @param marker_distribution A vector containing the position of all markers in Morgan.
#' @return The number of generations passed since the onset of hybridization
#' @examples
#' markers <- seq(from = 0, to = 1, length.out = 100)
#' J <- number_of_junctions_markers(N = 100, H_0 = 0.5, t = 200, marker_distribution = markers)
#' estimate_time_markers(J = J,
#'                      N = 100,
#' H_0 = 0.5,
#' marker_distribution = markers) #should be 200 again
#' @export
estimate_time_markers <- function(J = NA,
                                  N = Inf,
                                  H_0 = 0.5,
                                  marker_distribution = NA) {

 to_fit <- function(params) {
    expected_j <-
      number_of_junctions_markers(N = N,
                                  H_0 = H_0,
                                  t = params[[1]],
                                  marker_distribution)
    return(abs(expected_j - J))
  }

  upper_lim = 1e5
  if(J > 1e4) {
    upper_lim = J * 20
  }

  fitted <- stats::optimize(to_fit, interval = c(2, upper_lim) )
  return(fitted$minimum)
}