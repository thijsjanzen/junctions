#' Individual Based Simulation of the accumulation of junctions, given phasing
#' error
#' @description Individual based simulation of the accumulation of junctions,
#' returning phased and unphased data. Phased data is generated using a
#' user-provided error rate. Ancestry on both chromosomes of 10 randomly sampled
#' individuals per generations is returned.
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
#' recorded to be returned at the end of the simulation. If left at -1,
#' ancestry is recorded at every generation (computationally heavy).
#' @param seed Seed of the pseudo-random number generator
#' @param num_threads number of threads
#' @param verbose displays a progress bar
#' @param num_indiv_sampled the number of individuals sampled at each time point
#' to be genotyped
#' @param coverage fraction of markers that can be succesfully phased
#' @param error_rate fraction of markers that are erroneously
#' phased (e.g. swapped)
#' @return a tibble with five columns: [time, individual, marker location,
#'                                 ancestry chromosome 1, ancestry chromosome 2]
#' @examples
#' sim_phased_with_error(pop_size = 100, freq_ancestor_1 = 0.5,
#'                     total_runtime = 10, size_in_morgan = 1,
#'                     markers = 10, time_points = c(0, 5, 10),
#'                     seed = 42, coverage = 0.9, error_rate = 0.01)
#' @export
sim_phased_with_error <- function(pop_size = 100,
                                  freq_ancestor_1 = 0.5,
                                  total_runtime = 100,
                                  size_in_morgan = 1,
                                  markers = 100,
                                  time_points = -1,
                                  seed = NULL,
                                  num_threads = -1,
                                  verbose = TRUE,
                                  num_indiv_sampled = 10,
                                  coverage = 1,
                                  error_rate = 0) {

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

  markers <- get_num_markers(markers)

  record_true_junctions <- FALSE
  sim_output <-  sim_phased_unphased_cpp(pop_size,
                                         freq_ancestor_1,
                                         total_runtime,
                                         size_in_morgan,
                                         markers,
                                         time_points,
                                         seed,
                                         verbose,
                                         record_true_junctions,
                                         num_indiv_sampled,
                                         num_threads)

  true_data <- sim_output$results
  colnames(true_data) <- c("time", "individual", "location",
                           "anc_chrom_1", "anc_chrom_2")
  true_data <- tibble::as_tibble(true_data)

  markers <- unique(true_data$location)
  selected_markers <- sort(sample(markers,
                                  size = coverage * length(markers),
                                  replace = F))

  phased_data <- true_data[true_data$location %in% selected_markers, ]

  for (i in unique(phased_data$individual)) {
    focal_indices <- which(phased_data$individual == i)
    # select which ones to flip:

    sample_size <- stats::rbinom(size = length(focal_indices),
                                 n = 1,
                                 prob = error_rate)

    to_flip <- sample(focal_indices,
                      size = sample_size,
                      replace = F)

    if (length(to_flip) > 0) {
      temp <- phased_data$anc_chrom_1[to_flip]
      phased_data$anc_chrom_1[to_flip] <- phased_data$anc_chrom_2[to_flip]
      phased_data$anc_chrom_2[to_flip] <- temp
    }
  }

  return(list("true_data" = true_data,
              "phased_data" = phased_data))
}
