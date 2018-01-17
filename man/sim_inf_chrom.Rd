\name{sim_inf_chrom}
\alias{sim_inf_chrom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Individual Based Simulation of the accumulation of junctions
}
\description{
  Individual based simulation of the accumulation of junctions for a chromosome with an infinite number of recombination sites.
}
\usage{
sim_inf_chrom(pop_size, initial_heterozygosity, total_runtime, morgan, markers, seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pop_size}{
    Population Size
}
  \item{initial_heterozygosity}{
    Initial frequency of heterozygosity (e.g. in Hardy-Weinberg this equates to 2pq)
}
  \item{total_runtime}{
    Maximum time after which the simulation is to be stopped
}
  \item{morgan}{
    Size of the chromosome in Morgan
}
  \item{markers}{
    The number of genetic markers superimposed on the chromosome. If markers is set to -1, no markers are superimposed (faster simulation)
}
  \item{seed}{
    Seed of the pseudo-random number generator
}
}

\value{
  \item{avgJunctions}{vector of the average number of junctions at time = [0, total_runtime]}
  \item{detectedJunctions}{vector of the detected number of junctions at time = [0, total_runtime]}

}
\author{
Thijs Janzen
}

\examples{
  N = 100 #population size
  H_0 = 0.5 #initial heterozygosity
  maxT = 1000 #run time
  C = 1 #number of recombinations per meiosis
  numMarkers = 100
  seed = 42

  v <- sim_inf_chrom(N, H_0, maxT, C, numMarkers, seed)

  plot(v$avgJunctions, type = "l", xlab = "Generations",
        ylab = "Number of Junctions", main = "Example Infinite Chromosome")
  lines(v$detectedJunctions, col = "blue")
  legend("bottomright", c("Real number","Number detected"),
       lty = 1, col = c("black", "blue"))
}

\keyword{infinite chromosome}
\keyword{simulation}

