#' Individual Based Simulation of the accumulation of junctions
#' @description Individual based simulation of the accumulation of junctions, returning phased and unphased data. Ancestry on both chromosomes of 10 randomly sampled individuals per generations is returned.
#' @param pop_size Population Size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param total_runtime Maximum time after which the simulation is to be stopped
#' @param size_in_morgan Mean number of crossovers per meiosis (e.g. size in Morgan of the chromosome)
#' @param number_of_markers The number of genetic markers superimposed on the chromosome.
#' @param time_points vector with time points at which local ancestry has to be recorded to be returned at the end of the simulation. If left at -1, ancestry is recorded at every generation (computationally heavy).
#' @param seed Seed of the pseudo-random number generator
#' @param verbose displays a progress bar
#' @param num_threads if larger than one, multithreading is used.
#' @return a tibble with five columns: [time, individual, marker location, ancestry chromosome 1, ancestry chromosome 2]
#' @examples
#' sim_phased_unphased(pop_size = 100, freq_ancestor_1 = 0.5,
#'                     total_runtime = 10, size_in_morgan = 1,
#'                     number_of_markers = 10, time_points = c(0, 5, 10),
#'                     seed = 42)
#' @export
sim_phased_unphased <- function(pop_size = 100,
                                freq_ancestor_1 = 0.5,
                                total_runtime = 100,
                                size_in_morgan = 1,
                                number_of_markers = 100,
                                time_points = -1,
                                seed = NULL,
                                verbose = TRUE,
                                num_threads = 1) {
  if (length(time_points) == 1) {
    if (time_points == -1) {
      time_points <- seq(0, total_runtime, by = 1)
    }
  }

  if (is.null(seed)) {
    warning("warning! you did not provide a seed\n
            will use the time as a seed\n")
    seed <- Sys.time()
  }

  output <- sim_phased_unphased_cpp(pop_size,
                                      freq_ancestor_1,
                                      total_runtime,
                                      size_in_morgan,
                                      number_of_markers,
                                      time_points,
                                      seed,
                                      verbose,
                                      num_threads)

  colnames(output$results) <- c("time", "individual", "location",
                                "anc_chrom_1", "anc_chrom_2")

  return(tibble::as_tibble(output$results))
}
