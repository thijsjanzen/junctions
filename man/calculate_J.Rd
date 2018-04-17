\name{number_of_junctions}
\alias{number_of_junctions}

\title{
  Calculate the average number of junctions
}
\description{
  Calculate the average number of junctions after t generations, provided information on the initial heterozygosity, population size and the number of generations.
}
\usage{
  number_of_junctions(N = Inf, R = Inf, H_0 = 0.5, C = 1, t = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
    Population size
}
  \item{R}{
    Number of genetic markers
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
}

\value{
  The number of junctions for at time = t
}

\author{
Thijs Janzen
}

\examples{
  jt <-  number_of_junctions(N = 100, R = 1000, H_0 = 0.5, C = 1, t = 1000)
  jt2 <- number_of_junctions(N = 100, R = 1000, H_0 = 0.5, C = 1, t = 0:1000)
}

\keyword{junctions}
\keyword{analytic}