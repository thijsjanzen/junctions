\name{calculate_mat}
\alias{calculate_mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to calculate the maximum accurate time
}
\description{
  Function that calculates the maximum time after hybridization after which the number of junctions can still be reliably used to estimate the onset of hybridization. This is following equation 15 in Janzen et al. 2018.
}
\usage{
calculate_mat(N = Inf, R = Inf, H_0 = 0.5, C = 1)
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
}
\value{
  The maximum accurate time
}
\author{
  Thijs Janzen
}


\examples{
  calculate_mat(N = Inf, R = 1000, H_0 = 0.5, C = 1)
}

\keyword{analytic}
\keyword{maximum accurate time}