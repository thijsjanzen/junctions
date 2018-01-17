\name{calc_k}
\alias{calc_k}

\title{
  Calculate the limit of the number of junctions
}
\description{
  Calculate the average number of junctions after an infinite number of generations, provided information on the initial heterozygosity, population size and the number of generations.
}
\usage{
  calc_k(N = Inf, R = Inf, H_0 = 0.5, C = 1)
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
  The number of junctions for at time = infinity
}

\author{
  Thijs Janzen
}

\examples{
  k <-  calc_k(N = 100, R = 1000, H_0 = 0.5, C = 1)
}

\keyword{junctions}
\keyword{analytic}