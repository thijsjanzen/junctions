\name{estimate_time_markers}
\alias{estimate_time_markers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimate the time since the onset of hybridization, using the observed number of junctions, taking into account the distribution of markers
}
\description{
  Estimate the time since the onset of hybridization, following equation 1 in Janzen et al. unpublished
}
\usage{
estimate_time_markers(J = NA, N = Inf, H_0 = 0.5, marker_distribution = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{J}{
    The number of junctions
}
   \item{N}{
    Population size
}
   \item{H_0}{
    Initial frequency of heterozygosity (e.g. in Hardy-Weinberg this equates to 2pq)
}
  \item{marker_distribution}{
    A vector containing location of each marker in Morgan.
}

}
\value{
  The number of generations passed since the onset of hybridization
}
\author{
Thijs Janzen
}

\examples{
markers <- seq(from = 0, to = 1, length.out = 100)
J <- number_of_junctions_markers(N = 100, H_0 = 0.5, t = 200, marker_distribution = markers)
estimate_time_markers(J = J,
                      N = 100,
                      H_0 = 0.5,
                      marker_distribution = markers) #should be 200 again

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ time }% use one of  RShowDoc("KEYWORDS")
\keyword{ estimate }% __ONLY ONE__ keyword per line
\keyword{ analytic}
\keyword{ markers}