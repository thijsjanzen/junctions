#' Individual Based Simulation of the accumulation of junctions
#' @description Individual based simulation of the accumulation of junctions
#' for a chromosome with regularly distributed markers.
#' @param pop_size Population Size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param morgan Mean number of crossovers per meiosis (e.g. size in Morgan of
#' the chromosome)
#' @param total_runtime Maximum time after which the simulation is to be stopped
#' @param seed Seed of the pseudo-random number generator
#' @param R Number of regularly distributed markers
#' @return \item{avgJunctions}{vector of the average number of junctions at
#' time = [0, total_runtime]}
#' @examples
#' cat("example sim_fin_chrom")
#' sim_fin_chrom(pop_size = 100, freq_ancestor_1 = 0.5,
#'                    total_runtime = 10, morgan = 1, seed = 42,
#'                    R = 100)
#' @export
sim_fin_chrom <- function(pop_size = 100,
                           freq_ancestor_1 = 0.5,
                           total_runtime = 100,
                           morgan = 1,
                           seed = 42,
                           R = 100) {   # nolint
    set.seed(seed)
    .Call("_junctions_sim_fin_chrom", PACKAGE = "junctions",
        pop_size, freq_ancestor_1, total_runtime,
        morgan, seed, R)
}
