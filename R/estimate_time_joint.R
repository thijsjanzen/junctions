#' estimates the time since admixture, given unphased ancestry data.
#' @description Calculates the time since admixture, given unphased
#' ancestry data.
#' @param ancestry_information a matrix with five columns: column 1) indicator
#' of individual, column 2) indicator of chromosome, 3) location of marker in
#' Morgan, 4) ancestry at chromosome 5) ancestry at chromosome 2.
#' @param analysis_type how should the data be broken down? there are multiple
#' options: "all" - each chromosome is analysed separately, "individuals" - time
#' is jointly inferred for all chromosomes belonging to the same individual,
#' "chromosomes" - time is jointly inferred for all individuals with the same
#' chromosome.
#' @param phased is the data phased?
#' @param pop_size population size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param lower_lim lower limit of the optimization algorithm. Increase
#' if the expected admixture time is relatively ancient
#' @param upper_lim upper limit of hte optimization algorithm. If set too
#' large, recent admixture events can be overlooked - best to set as low
#' as possible.
#' @param num_threads num_threads, default is all threads. 5 threads is recommended.
#' @param verbose display intermediate output? Default = FALSE
#' @export
estimate_time_joint <- function(ancestry_information,
                                analysis_type = "all",
                                phased = FALSE,
                                pop_size = 1000,
                                freq_ancestor_1 = 0.5,
                                lower_lim = 2,
                                upper_lim = 2000,
                                num_threads = 1,
                                verbose = FALSE) {

  if (!(analysis_type %in% c("all", "individuals", "chromosomes"))) {
    stop("analysis type not known, did you perhaps spell individual instead of individuals?")
  }

  time_estimates <- c()
  if (analysis_type == "all") {
    for (indiv in unique(ancestry_information[, 1])) {
      for (chrom in unique(ancestry_information[, 2])) {
        local_anc_data <- subset(ancestry_information,
                                 ancestry_information[, 1] == indiv &
                                 ancestry_information[, 2] == chrom)

        result <- estimate_time_cpp(local_anc_matrix =
                                      as.matrix(local_anc_data[,c(2, 4, 5)]),
                                    locations = local_anc_data[, 3],
                                    pop_size = pop_size,
                                    freq_ancestor_1 = freq_ancestor_1,
                                    lower_lim = lower_lim,
                                    upper_lim = upper_lim,
                                    verbose = verbose,
                                    phased = phased,
                                    num_threads = num_threads)
        time_estimates <- rbind(time_estimates, c(indiv, chrom,
                                                  result$time, result$likelihood))
      }
    }
    colnames(time_estimates) <- c("individual", "chromosome",
                                  "time", "likelihood")
  }
  if (analysis_type == "individuals") {
    for (indiv in unique(ancestry_information[, 1])) {
      local_anc_data <- subset(ancestry_information,
                               ancestry_information[, 1] == indiv)

      result <- estimate_time_cpp(local_anc_matrix =
                                    as.matrix(local_anc_data[, c(2, 4, 5)]),
                                  locations = local_anc_data[, 3],
                                  pop_size = pop_size,
                                  freq_ancestor_1 = freq_ancestor_1,
                                  lower_lim = lower_lim,
                                  upper_lim = upper_lim,
                                  verbose = verbose,
                                  phased = phased,
                                  num_threads = num_threads)
      time_estimates <- rbind(time_estimates, c(indiv,
                                                result$time, result$likelihood))
    }
    colnames(time_estimates) <- c("individual", "time", "likelihood")
  }
  if (analysis_type == "chromosomes") {
    for (chrom in unique(ancestry_information[, 2])) {
      local_anc_data <- subset(ancestry_information,
                               ancestry_information[, 2] == chrom)

      result <- estimate_time_cpp(local_anc_matrix =
                                    as.matrix(local_anc_data[, c(1, 4, 5)]),
                                  locations = local_anc_data[, 3],
                                  pop_size = pop_size,
                                  freq_ancestor_1 = freq_ancestor_1,
                                  lower_lim = lower_lim,
                                  upper_lim = upper_lim,
                                  verbose = verbose,
                                  phased = phased,
                                  num_threads = num_threads)
      time_estimates <- rbind(time_estimates, c(chrom,
                                                result$time, result$likelihood))
    }
    colnames(time_estimates) <- c("chromosome", "time", "likelihood")
  }


  output <- tibble::as_tibble(time_estimates)
  return(output)
}