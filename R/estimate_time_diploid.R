#' estimates the time since admixture, given diploid ancestry data.
#' @description Calculates the time since admixture, given unphased
#' ancestry data.
#' @param ancestry_information a matrix with five columns: column 1) indicator
#' of individual, column 2) indicator of chromosome, 3) location of marker in
#' Morgan, 4) ancestry at chromosome 5) ancestry at chromosome 2.
#' @param analysis_type how should the data be broken down? there are multiple
#' options: "individuals" - time is inferred for each individual separately,
#' grouping all chromosomes together that belong to the same individual.
#' "chromosomes" - time is inferred for each chromosome separately, grouping
#' chromosomes together belonging from separate individuals. "separate" - time
#' is inferred for each chromosome from each individual separately, "all" - time
#' is inferred jointly for all chromosomes and individuals, grouping all
#' chromosomes and individuals together.
#' @param phased is the data phased?
#' @param pop_size population size
#' @param freq_ancestor_1 Frequency of ancestor 1 at t = 0
#' @param lower_lim lower limit of the optimization algorithm. Increase
#' if the expected admixture time is relatively ancient
#' @param upper_lim upper limit of hte optimization algorithm. If set too
#' large, recent admixture events can be overlooked - best to set as low
#' as possible.
#' @param num_threads num_threads, default is all threads. 5 threads is
#' recommended.
#' @param verbose display intermediate output? Default = FALSE
#' @return a tibble with the estimated time and its log likelihood. Depending
#' on analysis_type, columns identifying the individual and/or chromosome are
#' included, using the identifiers from the input. For analysis_type "all" only
#' time and loglikelihood are returned.
#' @export
estimate_time_diploid <- function(ancestry_information,
                                  analysis_type = "individuals",
                                  phased = FALSE,
                                  pop_size = 1000,
                                  freq_ancestor_1 = 0.5,
                                  lower_lim = 2,
                                  upper_lim = 2000,
                                  num_threads = 1,
                                  verbose = FALSE) {

  if (length(analysis_type) != 1 ||
      !(analysis_type %in% c("all", "individuals", "chromosomes", "separate"))) {
    stop("analysis_type must be one of \"all\", \"individuals\", ",
         "\"chromosomes\" or \"separate\"")
  }
  if (lower_lim >= upper_lim) {
    stop("lower_lim has to be smaller than upper_lim")
  }

  input <- prepare_diploid_input(ancestry_information)
  m <- input$data
  RcppParallel::setThreadOptions(num_threads)

  # group_col is passed to C++ as the chromosome identifier: C++ starts a new
  # chromosome whenever this value changes (or the location decreases).
  fit <- function(rows, group_col) {
    local_anc_data <- m[rows, , drop = FALSE]
    estimate_time_cpp(local_anc_matrix =
                        cbind(group_col[rows], local_anc_data[, 4:5,
                                                              drop = FALSE]),
                      locations = local_anc_data[, 3],
                      pop_size = pop_size,
                      freq_ancestor_1 = freq_ancestor_1,
                      lower_lim = lower_lim,
                      upper_lim = upper_lim,
                      verbose = verbose,
                      phased = phased)
  }

  indiv_code <- m[, 1]
  chrom_code <- m[, 2]
  # unique id for each chromosome of each individual, see
  # prepare_diploid_input
  chrom_id <- m[, 6]
  # id for each individual / chromosome label combination
  label_id <- match(paste(indiv_code, chrom_code),
                    unique(paste(indiv_code, chrom_code)))

  if (analysis_type == "all") {
    result <- fit(seq_len(nrow(m)), chrom_id)
    return(tibble::tibble(time = result$time,
                          loglikelihood = result$likelihood))
  }

  if (analysis_type == "individuals") {
    groups <- unique(indiv_code)
    results <- lapply(groups, function(i) fit(which(indiv_code == i),
                                              chrom_id))
    output <- tibble::tibble(individual = input$indiv_levels[groups])
  }

  if (analysis_type == "chromosomes") {
    groups <- sort(unique(chrom_code))
    results <- lapply(groups, function(i) fit(which(chrom_code == i),
                                              chrom_id))
    output <- tibble::tibble(chromosome = input$chrom_levels[groups])
  }

  if (analysis_type == "separate") {
    groups <- unique(label_id)
    results <- lapply(groups, function(i) fit(which(label_id == i),
                                              chrom_id))
    first_row <- match(groups, label_id)
    output <- tibble::tibble(
      individual = input$indiv_levels[indiv_code[first_row]],
      chromosome = input$chrom_levels[chrom_code[first_row]])
  }

  output$time <- vapply(results, function(r) r$time, numeric(1))
  output$loglikelihood <- vapply(results, function(r) r$likelihood,
                                 numeric(1))
  return(output)
}
