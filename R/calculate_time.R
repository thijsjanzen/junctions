estimate_time <- function(J = NA, 
                          N = Inf, 
                          R = Inf, 
                          H_0 = 0.5, 
                          C = 1) {
  if (is.na(J)) {
    cat("ERROR! did you forget to provide J?")
    return()
  }
  
  if (is.infinite(N) && is.infinite(R)) {
    cat("both N and R are infinite\n")
    cat("can not estimate t\n")
  }
  
  K <- calc_k(N, R, H_0, C)

  u <- 1 - 1 / (2 * N) - C / R
  
  t <- log(1 - J / K) / (log(u))
  return(t)
}
