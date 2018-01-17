\name{sim_fin_chrom}
\alias{sim_fin_chrom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Individual Based Simulation of the accumulation of junctions
}
\description{
  Individual based simulation of the accumulation of junctions for a chromosome with an finite number of recombination sites.
}
\usage{
sim_fin_chrom(pop_size, initial_heterozygosity, total_runtime, morgan, seed, R)
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

  \item{seed}{
    Seed of the pseudo-random number generator
  }

  \item{R}{
    Number of genetic markers
  }
}

\value{
  \item{avgJunctions}{vector of the average number of junctions at time = [0, total_runtime]}
}
\author{
Thijs Janzen
}

\examples{
  N = 100 #population size
  R = 100
  H_0 = 0.5 #initial heterozygosity
  maxT = 1000 #run time
  C = 1 #number of recombinations per meiosis
  seed = 42

  v <- sim_fin_chrom(N, H_0, maxT, C, seed, R)
}

\keyword{finite chromosome}
\keyword{simulation}