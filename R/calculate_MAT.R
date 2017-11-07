calculate_MAT <-
function(N = Inf, R = Inf, H_0 = 0.5, C = 1) {
  if(is.infinite(N) && is.infinite(R)) {
    cat("both N and R are infinite\n")
    cat("can not estimate MAT\n")
  }

  K <- calc_K(N, R, H_0, C)
  
  u <- 1 - 1/(2*N) - C/R
  MAT = log(1/K)/log(u)
  return(MAT)
}
