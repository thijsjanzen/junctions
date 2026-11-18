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
  args <- recycle_args(N, R, H_0, C)
  N <- args[[1]]  # nolint
  R <- args[[2]]  # nolint
  H_0 <- args[[3]]  # nolint
  C <- args[[4]]  # nolint

  # variable names are consistent with Janzen et al. 2018, thus nolint
  K <- H_0 * C * 2 * N * R / (2 * N * C + R)  # nolint

  inf_n <- is.infinite(N)
  inf_r <- is.infinite(R)
  K[inf_n] <- H_0[inf_n] * R[inf_n]                    # nolint
  K[inf_r] <- H_0[inf_r] * C[inf_r] * 2 * N[inf_r]     # nolint
  K[inf_n & inf_r] <- Inf                              # nolint
  return(K)
}
