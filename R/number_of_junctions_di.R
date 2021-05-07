#' Calculate the expected number of junctions between two markers separated by
#' a given amount of recombination
#' @description Calculate the expected number of junctions after t generations,
#'  provided information on the initial heterozygosity, population size, the
#'  number of generations since the onset of admixture and the distance between
#'  two markers.
#' @param N Population Size
#' @param H_0 Frequency of heterozygosity at t = 0
#' @param t Time since admixture
#' @param di Distance between two markers in Morgan
#' @return Estimated number of junctions at time t
#' @examples
#' number_of_junctions_di(N = 100, H_0 = 0.5, t = 1000, di = 0.01)+
#' @export
number_of_junctions_di <- function(N = Inf,        # nolint
                                   H_0 = 0.5,      # nolint
                                   t = 100,
                                   di = 1e-6) {
  local_k <- 2 * N  * H_0 / (2 * N  + 1 / di)
  if (is.infinite(N)) local_k <- H_0

  local_j <- local_k - local_k * (1 - 1 / (2 * N) - di) ^ t
  return(local_j)
}
