calc_k <- function(N = Inf, 
                   R = Inf, 
                   H_0 = 0.5, 
                   C = 1) {
  
  K <- H_0 * C * 2 * N * R / (2 * N * C + R);

  if (is.infinite(N)) {
    K <- H_0 * R
  }

  if (is.infinite(R)) {
    K <- H_0 * C * 2 * N
  } 

  if (is.infinite(N) && is.infinite(R)) {
    K <- Inf
  }
  return(K)
}