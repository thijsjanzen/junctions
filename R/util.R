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
    return(markers)
  }
}
# nolint end