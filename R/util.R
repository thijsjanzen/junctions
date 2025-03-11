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

    if (max(markers) != 1) {
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

    to_flip <- sample(focal_indices,
                      size = sample_size,
                      replace = FALSE)

    if (length(to_flip) > 0) {
      temp <- phased_data$anc_chrom_1[to_flip]
      phased_data$anc_chrom_1[to_flip] <- phased_data$anc_chrom_2[to_flip]
      phased_data$anc_chrom_2[to_flip] <- temp
    }
  }

  return(list("true_data" = true_data,  "phased_data" = phased_data))
}

check_time_points <- function(time_points, max_time) {

  within_range <- time_points <= max_time
  time_points <- time_points[within_range]
  if (length(time_points) < 1) {
    warning("all chosen time points were past the simulation time,
             chosen to only measure at the last time step")
    time_points <- max_time
  }

  if (length(time_points) == 1) {
    if (time_points[[1]] == -1) {
      time_points <- seq(0, max_time, by = 1)
    }
  }
  return(time_points)
}

# nolint end
