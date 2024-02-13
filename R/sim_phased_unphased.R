#' Individual Based Simulation of the accumulation of junctions
#' @description Individual based simulation of the accumulation of junctions,
#' returning phased and unphased data. Ancestry on both chromosomes of 10
#' randomly sampled individuals per generations is returned.
#' @param pop_size Population Size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param total_runtime Maximum time after which the simulation is to be stopped
#' @param size_in_morgan Mean number of crossovers per meiosis (e.g. size in
#' Morgan of the chromosome)
#' @param markers If a single number is provided, the number is used as the
#' total number of markers generated either randomly, or using a regular
#' distribution (a regular distribution is chosen if the number is negative). If
#' a vector is provided, that vector is used.
#' @param time_points vector with time points at which local ancestry has to be
#'  recorded to be returned at the end of the simulation. If left at -1,
#'  ancestry is recorded at every generation (computationally heavy).
#' @param num_threads default is 1. -1 takes all available threads.
#' @param verbose displays a progress bar
#' @param record_true_junctions keep track of the true number of junctions?
#' @param num_indiv_sampled the number of individuals sampled at each time point
#' to be genotyped
#' @param coverage fraction of markers that can be succesfully phased
#' @param error_rate fraction of markers that are erroneously
#' phased (e.g. swapped)
#' @return a tibble with five columns: [time, individual, marker location,
#'                             ancestry chromosome 1, ancestry chromosome 2]
#' @examples
#' \dontrun{
#' sim_phased_unphased(pop_size = 100, freq_ancestor_1 = 0.5,
#'                     total_runtime = 10, size_in_morgan = 1,
#'                     markers = 10, time_points = c(0, 5, 10),
#'                     num_threads = 1)
#' }
#' @export
sim_phased_unphased <- function(pop_size = 100,
                                freq_ancestor_1 = 0.5,
                                total_runtime = 100,
                                size_in_morgan = 1,
                                markers = 100,
                                time_points = -1,
                                num_threads = 1,
                                verbose = FALSE,
                                record_true_junctions = FALSE,
                                num_indiv_sampled = 10,
                                coverage = 1,
                                error_rate = 0) {

  RcppParallel::setThreadOptions(num_threads)

  time_points <- check_time_points(time_points, total_runtime)

  markers <- get_num_markers(markers)

  output <- sim_phased_unphased_cpp(pop_size,
                                    freq_ancestor_1,
                                    total_runtime,
                                    size_in_morgan,
                                    markers,
                                    time_points,
                                    verbose,
                                    record_true_junctions,
                                    num_indiv_sampled,
                                    num_threads)

  if (coverage != 1 || error_rate != 0) {
    phasing_result <- apply_phasing_error(output, coverage, error_rate)

    return(list("true_data" = phasing_result$true_data,
                "phased_data" = phasing_result$phased_data))
  } else {
    colnames(output$results) <- c("time", "individual", "location",
                                  "anc_chrom_1", "anc_chrom_2")

    if (!record_true_junctions)
      return(tibble::as_tibble(output$results))

    if (record_true_junctions) {
      colnames(output$true_results) <- c("time",
                                         "individual",
                                         "junctions_chrom_1",
                                         "junctions_chrom_2")
      output <- list("results" = tibble::as_tibble(output$results),
                     "true_results" = tibble::as_tibble(output$true_results))
      return(output)
    }
  }
}
