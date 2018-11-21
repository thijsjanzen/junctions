\name{number_of_junctions_markers}
\alias{number_of_junctions_markers}

\title{
  Calculate the expected number of junctions given a distribution of markers
}
\description{
  Calculate the expected number of junctions after t generations, provided information on the initial heterozygosity, population size, the number of generations since the onset of admixture and the distribution of markers.
}
\usage{
  number_of_junctions_markers(N = Inf, H_0 = 0.5, C = 1, t = 100, marker_distribution = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
    Population size
}
  \item{H_0}{
    Initial frequency of heterozygosity (e.g. in Hardy-Weinberg this equates to 2pq)
}
  \item{C}{
    Mean number of crossovers per meiosis (e.g. size in Morgan of the chromosome)
}
  \item{t}{
    A single time point or a vector of time points for which the number of junctions is estimated.
}
  \item{marker_distribution}{
    A vector containing the relative position (e.g. in [0, 1]) of all markers. Please note that this is NOT the distribution in Morgans!
}
}

\value{
  The number of junctions for at time = t
}

\author{
Thijs Janzen
}

\examples{
  markers <- seq(from = 0, to = 1, length.out = 1000)
  jt <-  number_of_junctions_markers(N = 100, 
                                     H_0 = 0.5, 
                                     C = 1, 
                                     t = 1000, 
                                     marker_distribution = markers)
  random_markers <- sort(runif(1000, 0, 1))
  jt2 <- number_of_junctions_markers(N = 100, 
                                     H_0 = 0.5, 
                                     C = 1, 
                                     t = 1000, 
                                     marker_distribution = random_markers)
}

\keyword{junctions}
\keyword{analytic}