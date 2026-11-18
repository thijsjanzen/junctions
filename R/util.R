# nolint start
#' @keywords internal
#' @rawNamespace useDynLib(junctions)
#' @rawNamespace import(Rcpp)
#' @rawNamespace import(nloptr)
#' @rawNamespace importFrom(RcppParallel, RcppParallelLibs)
get_num_markers <- function(markers) {
  if (length(markers) == 1) {
    num_markers <- abs(markers[1])

    if (markers[1] < 0) { # evenly spaced markers
      di = 1.0 / (num_markers);
      markers <- seq(di, 1 - di, length.out = num_markers)
      return(markers)
    } else {
      markers <- sort(stats::runif(n = num_markers, min = 0, max = 1))
      markers <- unique(markers)
      return(markers)
    }
  } else {
    if (anyNA(markers) || any(markers < 0)) {
      stop("marker positions must be non-negative and not NA")
    }
    markers <- sort(unique(markers))
    # positions are relative to chromosome length; only rescale when they
    # exceed 1, so user-supplied positions in [0, 1] are used as given
    if (max(markers) > 1) {
      markers <- markers / max(markers)
    }

    return(markers)
  }
}

#' @keywords internal
apply_phasing_error <- function(output,
                                coverage,
                                error_rate) {

  true_data <- output$results
  colnames(true_data) <- c("time", "individual", "location",
                           "anc_chrom_1", "anc_chrom_2")
  true_data <- tibble::as_tibble(true_data)

  markers <- unique(true_data$location)
  selected_markers <- sort(sample(markers,
                                  size = coverage * length(markers),
                                  replace = FALSE))

  phased_data <- true_data[true_data$location %in% selected_markers, ]

  for (i in unique(phased_data$individual)) {
    focal_indices <- which(phased_data$individual == i)
    # select which ones to flip:
    sample_size <- stats::rbinom(size = length(focal_indices),
                                 n = 1,
                                 prob = error_rate)

    to_flip <- focal_indices[sample.int(length(focal_indices), sample_size)]

    if (length(to_flip) > 0) {
      temp <- phased_data$anc_chrom_1[to_flip]
      phased_data$anc_chrom_1[to_flip] <- phased_data$anc_chrom_2[to_flip]
      phased_data$anc_chrom_2[to_flip] <- temp
    }
  }

  return(list("true_data" = true_data,  "phased_data" = phased_data))
}

check_time_points <- function(time_points, max_time) {

  if (length(time_points) == 1 && time_points == -1) {
    return(seq(0, max_time, by = 1))
  }
  if (any(time_points < 0 | time_points != round(time_points))) {
    stop("time_points must be non-negative whole numbers (or -1 for all)")
  }

  within_range <- time_points <= max_time
  time_points <- time_points[within_range]
  if (length(time_points) < 1) {
    warning("all chosen time points were past the simulation time, ",
            "chosen to only measure at the last time step")
    time_points <- max_time
  }

  return(time_points)
}

# evaluate expr with R's random number generator seeded by `seed`, restoring
# the user's RNG state afterwards (so the global state is not affected)
with_local_seed <- function(seed, expr) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  expr
}

# recycle all arguments to a common length, so the analytical
# functions can be called with vectors for any of their arguments
recycle_args <- function(...) {
  args <- list(...)
  n <- max(lengths(args))
  lapply(args, rep_len, length.out = n)
}

# Checks and standardises diploid ancestry input. Individual and chromosome
# identifiers can be of any type (numeric, character, factor); they are
# converted to integer codes for the C++ code, and the original values are
# returned so the output can be labelled with them.
# Rows are sorted by individual, chromosome and position (keeping the order in
# which individuals and chromosomes first appear), so the C++ code sees each
# chromosome as one contiguous, increasing block.
prepare_diploid_input <- function(ancestry_information) {
  x <- as.data.frame(ancestry_information)
  if (ncol(x) != 5) {
    stop("ancestry_information must have five columns: individual, ",
         "chromosome, location, ancestry chromosome 1, ancestry chromosome 2")
  }

  indiv <- x[[1]]
  chrom <- x[[2]]
  pos <- as.numeric(x[[3]])
  anc <- cbind(as.numeric(x[[4]]), as.numeric(x[[5]]))

  if (anyNA(indiv) || anyNA(chrom) || anyNA(pos)) {
    stop("individual, chromosome and location columns can not contain NA")
  }
  check_ancestry_values(anc)

  indiv_levels <- unique(indiv)
  chrom_levels <- unique(chrom)
  indiv_code <- match(indiv, indiv_levels)
  chrom_code <- match(chrom, chrom_levels)

  m <- cbind(indiv_code, chrom_code, pos, anc[, 1], anc[, 2])
  m <- m[order(indiv_code, chrom_code, pos), , drop = FALSE]

  # chromosomes with a single marker contain no information on junctions
  group <- paste(m[, 1], m[, 2])
  markers_per_chrom <- table(group)
  single <- names(markers_per_chrom)[markers_per_chrom < 2]
  if (length(single) == length(markers_per_chrom)) {
    stop("each chromosome needs at least two markers to estimate time")
  }
  if (length(single) > 0) {
    warning(length(single), " chromosome(s) with a single marker ",
            "carry no information on junctions and were ignored")
    m <- m[!(group %in% single), , drop = FALSE]
  }

  colnames(m) <- NULL
  list(data = m,
       indiv_levels = indiv_levels,
       chrom_levels = chrom_levels)
}

check_ancestry_values <- function(anc) {
  if (anyNA(anc)) {
    stop("ancestry columns can not contain NA")
  }
  if (!all(anc %in% c(0, 1))) {
    stop("ancestry has to be coded as 0 (ancestor 1) and 1 (ancestor 2), ",
         "found values: ", paste(utils::head(setdiff(unique(c(anc)), c(0, 1))),
                                 collapse = ", "))
  }
}

# nolint end
