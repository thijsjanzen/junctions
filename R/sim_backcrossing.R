#' Function to simulate data using a back crossing scheme
#' @description Individual based simulation of the accumulation of junctions,
#' under a back crossing scheme
#' @param population_size Population size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param total_runtime Number of generations to simulate
#' @param size_in_morgan Mean number of crossovers per meiosis (e.g. size in
#' Morgan of the chromosome)
#' @param number_of_markers number of molecular markers
#' @param seed Seed of the pseudo-random number generator
#' @param time_points vector with time points at which local ancestry has to be
#' recorded to be returned at the end of the simulation. If left at -1,
#' ancestry is recorded at every generation (computationally heavy).
#' @return List with five entries: average_junctions: average number of
#' junctions over time, detected_junctions: average number of detected
#' junctions, given the markers. markers: vector with the locations of the
#' molecular markers, junction_distribution: distribution of junctions per
#' time step average_heterozygosity: average heterozygosity.
#' @examples
#' sim_backcrossing(population_size = 100,
#'                        total_runtime = 5,
#'                        size_in_morgan = 1,
#'                        number_of_markers = 100,
#'                        seed = 6,
#'                        time_points = 1:5)
#' @export
sim_backcrossing <- function(population_size = 100,
                             freq_ancestor_1 = 0.5,
                             total_runtime = 5,
                             size_in_morgan = 1,
                             number_of_markers = 100,
                             seed = 6,
                             time_points = -1) {
  if (length(time_points) == 1) {
    if (time_points == -1) {
      time_points <- seq(0, total_runtime, by = 1)
    }
  }

  output <- simulate_backcrossing_cpp(population_size,
                                      freq_ancestor_1,
                                      total_runtime,
                                      size_in_morgan,
                                      number_of_markers,
                                      time_points,
                                      seed)
  return(output)
}
