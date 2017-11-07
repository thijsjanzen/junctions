calculate_J <-
function(N = Inf, R = Inf, H_0 = 0.5, C = 1, t = 100) {
  ######## N: population size
  ######## R: number of evenly spaced recombination sites
  ######## H_0: initial heterozygosity 2pq (at t = 0)
  ######## C: size of the chromosome in Morgan
  ######## maxT: maximum number of timesteps
  ########
  ######## returns: number of junctions in t = [0, maxT], including ZERO!
  
  if(is.infinite(N) && is.infinite(R)) {
    # If both N and R are infinite, R gives 
    # numerical problems using equation 12 
    # to calculate K, so instead we use
    # equation 1
    jt <- H_0 * C * t;
    return(jt)
  }
  
  K <- calc_K(N, R, H_0, C)
  
  jt <- K - K *(1-H_0*C/K)^t
  return(jt) #jt in [0,maxT]
}
