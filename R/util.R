# nolint start
#' @keywords internal
#' @rawNamespace useDynLib(junctions)
#' @rawNamespace import(Rcpp)
#' @rawNamespace import(nloptr)
get_num_markers <- function(markers) {
  if (length(markers) == 1) {
    num_markers <- abs(markers[1])

    if (markers[1] < 0) { # evenly spaced markers
      di = 1.0 / (num_markers);
      markers <- seq(di, 1 - di, length.out = num_markers)
      return(markers)
    } else {
      markers <- stats::runif(n = num_markers)
      markers <- unique(markers)
      while (length(markers) < num_markers) {
        add <- stats::runif(num_markers - length(markers))
        markers <- unique(c(markers, add))
      }
      markers <- sort(markers)
      return(markers)
    }
  } else {
    return(markers)
  }
}

#' @keywords internal
single_state <- function(t, N, d) {

  trans_matrix <- matrix(0, 7, 7)
  trans_matrix[1, ] <- c(1 - 1 / (2*N) - 2 * d , 2 * d, 0, 0, 0, 1 / (2*N), 0)
  trans_matrix[2, ] <- c(1 / (2*N), 1 - 3 * 1 / (2*N) - d, d, 2 * 1 / (2*N), 0, 0, 0)
  trans_matrix[3, ] <- c(0, 2 * 1 / (2*N), 1 - 4 * 1/(2*N), 0, 2 * 1/(2*N), 0, 0)
  trans_matrix[4, ] <- c(0, 0, 0, 1 - 1 / (2*N) - d, d, 1/(2*N), 0)
  trans_matrix[5, ] <- c(0, 0, 0, 2 * 1 / (2*N), 1 - 3 * 1/(2*N), 0, 1 / (2*N))
  trans_matrix[6, ] <- c(0, 0, 0, 0, 0, 1 - d, d)
  trans_matrix[7, ] <- c(0 ,0, 0, 0, 0, 1 / (2*N),  1 - 1 / (2*N))

  initial_state <- c(1, 0, 0, 0, 0, 0, 0)

  output_state <- initial_state %*% expm::`%^%`(trans_matrix, t)
  return(output_state)
}

# nolint end