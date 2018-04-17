\name{estimate_time}
\alias{estimate_time}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimate the time since the onset of hybridization, using the number of junctions
}
\description{
  Estimate the time since the onset of hybridization, following equation 14 in Janzen et al. 2018
}
\usage{
estimate_time(J = NA, N = Inf, R = Inf, H_0 = 0.5, C = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{J}{
    The number of junctions
}
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
  The number of generations passed since the onset of hybridization
}
\author{
Thijs Janzen
}

\examples{
J <- number_of_junctions(N = 100, R = 1000, H_0 = 0.5, C = 1, t = 200)
estimate_time(J = J, N = 100, R = 1000, H_0 = 0.5, C = 1) #should be 200 again

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ time }% use one of  RShowDoc("KEYWORDS")
\keyword{ estimate }% __ONLY ONE__ keyword per line
\keyword{ analytic}
