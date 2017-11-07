\name{calculate_J}
\alias{calculate_J}

\title{
  Calculate the average number of junctions
}
\description{
  Calculate the average number of junctions after t generations, provided information on the initial heterozygosity, population size and the number of generations.
}
\usage{
  calculate_J(N = Inf, R = Inf, H_0 = 0.5, C = 1, t = 100)
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
    The number of junctions will be calculated for this timepoint.
}
}

\value{
  The number of junctions for at time = t
}
\references{
 Janzen, Thijs, Arne W. Nolte, and Arne Traulsen. "Haplotype block dynamics in hybrid populations." bioRxiv (2016): 058107.
}

\author{
Thijs Janzen
}

\examples{
  calculate_J(N = 100, R = 1000, H_0 = 0.5, C = 1, t = 1000)
  calculate_J(N = 100, R = 1000, H_0 = 0.5, C = 1, t = 0:1000)
}

