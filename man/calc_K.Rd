\name{calc_K}
\alias{calc_K}

\title{
  Calculate K
}
\description{
  Calculates the number of junctions at t = infinity, following equation 12 in Janzen et al. 2017.
}
\usage{
calc_K(N = Inf, R = Inf, H_0 = 0.5, C = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
        Population Size
}
  \item{R}{
        Number of Genetic Markers
}
  \item{H_0}{
        Initial frequency of heterozygosity (e.g. in HW: 2pq)
}
  \item{C}{
        Number of crossovers per chromosome (e.g. the size of the chromosome in Morgan)
}
}

\value{
  The value of K
}
\references{
 Janzen, Thijs, Arne W. Nolte, and Arne Traulsen. "Haplotype block dynamics in hybrid populations." bioRxiv (2016): 058107.
}
\author{
  Thijs Janzen
}

\examples{
  calc_K(N = 100, R = 100, H_0 = 0.5, C = 1)
}
