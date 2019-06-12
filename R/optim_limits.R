#' Optimizes the lower and upper limit for estimation of time
#' @description Optimizes the lower and upper limit for estimation of time
#' @param lower initial lower limit
#' @param upper initial upper limit
#' @param calc_func function to calculate the time, should accept two arguments: upper, lower and should return a stats::optimize object with two properties: $minimum and $objective
#' @param iterations number of times the boundaries should be adjusted to obtain an estimate that is not close to the boundaries. After this number of iterations is exceeded, NA is returned
#' @return final estimate
#' @export
optim_limits <- function(lower = 1,
                         upper = 1e4,
                         calc_func,
                         iterations = 10) {

  estim <- calc_func(lower, upper)

  is_close_to_lower <- lower / estim$minimum > 0.9
  is_close_to_upper <- estim$minimum / upper > 0.9

  lower_cnt <- 1
  cat(lower, upper, estim$minimum, "\n")
  while (is_close_to_lower) {

    lower <- max(0.5 * lower, 1)
    upper <- upper * 0.5
    estim <- calc_func(lower, upper)
    is_close_to_lower <- lower / estim$minimum > 0.9
    lower_cnt <- lower_cnt + 1
    cat(lower, upper, estim$minimum, "\n")
    if (lower_cnt >= iterations) break
  }

  upper_cnt <- 1
  while (is_close_to_upper) {
    upper <- upper * 1.5
    lower <- lower * 1.5
    estim <- calc_func(lower, upper)
    is_close_to_upper <- estim$minimum / upper > 0.9
    upper_cnt <- lower_cnt + 1
    cat(lower, upper, estim$minimum, "\n")
    if (upper_cnt >= iterations) break
  }

  is_close_to_lower <- lower / estim$minimum > 0.9
  is_close_to_upper <- estim$minimum / upper > 0.9

  if (is_close_to_lower || is_close_to_upper) {
    estim$minimum <- NA
    estim$objective <- NA
  }
  return(estim)
}