#' Calculate the average number of junctions
#' @description Calculate the average number of junctions in a single chromosome
#' after t generations, provided information on the initial heterozygosity,
#' population size and the number of generations.
#' @param N Population Size
#' @param R Number of genetic markers
#' @param H_0 Frequency of heterozygosity at t = 0
#' @param C Mean number of crossovers per meiosis (e.g. size in Morgan of the
#' chromosome)
#' @param t Time since admixture
#' @return Estimated number of junctions at time t
#' @examples
#' jt <-  number_of_junctions(N = 100, R = 1000, H_0 = 0.5, C = 1, t = 1000)
#' jt2 <- number_of_junctions(N = 100, R = 1000, H_0 = 0.5, C = 1, t = 0:1000)
#' @export
number_of_junctions <- function(N = Inf,      # nolint
                                R = Inf,      # nolint
                                H_0 = 0.5,    # nolint
                                C = 1,        # nolint
                                t = 100) {

  if (is.infinite(N) && is.infinite(R)) {
    # If both N and R are infinite, R gives
    # numerical problems using equation 12
    # to calculate K, so instead we use
    # equation 1
    jt <- H_0 * C * t
    return(jt)
  }

  K <- junctions::calc_k(N, R, H_0, C)   # nolint

  jt <- K - K * (1 - H_0 * C / K) ^ t
  return(jt)
}
