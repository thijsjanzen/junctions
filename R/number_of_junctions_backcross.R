#' Calculate the average number of junctions during backcrossing
#' @description Calculate the expected number of junctions after t
#' generations, in a backcrossing mating scheme.
#' @param H_0 Frequency of heterozygosity at t = 0
#' @param C Mean number of crossovers per meiosis (e.g. size in Morgan of the
#' chromosome)
#' @param t Time since admixture
#' @return Estimated number of junctions at time t
#' @examples
#' jt <-  number_of_junctions_backcross(H_0 = 0.1, C = 1, t = 5)
#' @export
number_of_junctions_backcross <- function(H_0 = 0.5,    # nolint
                                          C = 1,        # nolint
                                          t = 100) {
  jt <- 2 * H_0 * t * C * 2 ^ (-t)
  return(jt)
}
