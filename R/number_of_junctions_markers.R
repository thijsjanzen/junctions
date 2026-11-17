#' Calculate the expected total number of junctions in a chromosome, given the
#' distribution of markers
#' @description Calculate the expected number of junctions after t generations,
#' provided information on the initial heterozygosity, population size, the
#' number of generations since the onset of admixture and the distribution of
#' markers.
#' @param N Population Size
#' @param H_0 Frequency of heterozygosity at t = 0
#' @param t Time since admixture
#' @param marker_distribution A vector containing the position of all markers
#' in Morgan.
#' @return Estimated number of observed junctions at time t
#' @examples
#' markers <- seq(from = 0, to = 1, length.out = 1000)
#' jt <-  number_of_junctions_markers(N = 100,
#'                                   H_0 = 0.5,
#'                                   t = 1000,
#'                                   marker_distribution = markers)
#' random_markers <- sort(runif(1000, 0, 1))
#' jt2 <- number_of_junctions_markers(N = 100,
#'                                   H_0 = 0.5,
#'                                   t = 1000,
#'                                   marker_distribution = random_markers)
#' @export
number_of_junctions_markers <- function(N = Inf,       # nolint
                                        H_0 = 0.5,     # nolint
                                        t = 100,
                                        marker_distribution = NA) {

  if (length(marker_distribution) < 2) {
    stop("not provided a vector containing sufficient marker locations")
  }

  number_of_junctions <- 0
  for (i in 2:length(marker_distribution)) {
    di <- marker_distribution[i] - marker_distribution[i - 1]
    expected_junctions <- number_of_junctions_di(N, H_0, t, di)
    number_of_junctions <- number_of_junctions + expected_junctions
  }
  return(number_of_junctions)
}
