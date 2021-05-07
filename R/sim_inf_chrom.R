#' Individual Based Simulation of the accumulation of junctions
#' @description Individual based simulation of the accumulation of junctions
#' for a chromosome with an infinite number of recombination sites.
#' @param pop_size Population Size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param total_runtime Maximum time after which the simulation is to be stopped
#' @param morgan Mean number of crossovers per meiosis (e.g. size in Morgan of
#' the chromosome)
#' @param markers The number of genetic markers superimposed on the chromosome.
#' If markers is set to -1, no markers are superimposed (faster simulation)
#' @param seed Seed of the pseudo-random number generator
#' @return \item{avgJunctions}{vector of the average number of junctions at
#' time = [0, total_runtime]}
#' @examples
#' cat("example sim inf chrom")
#' v <- sim_inf_chrom(pop_size = 100, freq_ancestor_1 = 0.5,
#'                    total_runtime = 10, morgan = 1, markers = 100,
#'                    seed = 42)
#' plot(v$avgJunctions, type = "l", xlab = "Generations",
#' ylab = "Number of Junctions", main = "Example Infinite Chromosome")
#' lines(v$detectedJunctions, col = "blue")
#' legend("bottomright", c("Real number","Number detected"),
#'        lty = 1, col = c("black", "blue"))
#' @export
sim_inf_chrom <- function(pop_size = 100,
                          freq_ancestor_1 = 0.5,
                          total_runtime = 100,
                          morgan = 1,
                          markers = -1,
                          seed = 42) {
    set.seed(seed)
    .Call("_junctions_sim_inf_chrom", PACKAGE = "junctions",
        pop_size, freq_ancestor_1, total_runtime,
        morgan, markers, seed)
}
