#' Calculate the limit of the number of junctions
#' @description Calculate the average number of junctions after an infinite
#' number of generations, provided information on the initial heterozygosity,
#' population size and the number of generations.
#' @param N population size
#' @param R number of markers
#' @param H_0 initial heterozygosity (at the time of admixture)
#' @param C Mean number of crossovers per meiosis (e.g. size in Morgan of
#' the chromosome)
#' @return The number of junctions for at time = infinity
#' @examples
#' k <-  calc_k(N = 100, R = 1000, H_0 = 0.5, C = 1)
#' @keywords junctions
#' @export
calc_k <- function(N = Inf,   # nolint
                   R = Inf,   # nolint
                   H_0 = 0.5, # nolint
                   C = 1) {   # nolint

  # variable names are consistent with Janzen et al. 2018, thus nolint
  K <- H_0 * C * 2 * N * R / (2 * N * C + R)  # nolint

  if (is.infinite(N)) {
    K <- H_0 * R      # nolint
  }

  if (is.infinite(R)) {
    K <- H_0 * C * 2 * N   # nolint
  }

  if (is.infinite(N) && is.infinite(R)) {
    K <- Inf    # nolint
  }
  return(K)
}
