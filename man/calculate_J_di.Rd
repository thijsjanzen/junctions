\name{number_of_junctions_di}
\alias{number_of_junctions_di}

\title{
  Calculate the expected number of junctions between two markers separated by di
}
\description{
  Calculate the expected number of junctions after t generations, provided information on the initial heterozygosity, population size, the number of generations since the onset of admixture and the distance between two markers.
}
\usage{
  number_of_junctions_di(N = Inf, H_0 = 0.5, C = 1, t = 100, di = 1e-6)
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
  \item{di}{
    Distance between two markers, relative to the total length of the chromosome
}
}

\value{
  The number of junctions for at time = t
}

\author{
Thijs Janzen
}

\examples{
  number_of_junctions_di(N = 100, H_0 = 0.5, C = 1, t = 1000, di = 0.01)
}

\keyword{junctions}
\keyword{analytic}