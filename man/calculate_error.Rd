\name{time_error}
\alias{time_error}
\title{
  Estimate the error in the time estimate
}
\description{
  Calculate the error in the estimate of the onset of hybridisation, following Equations 3 & 4 in the Supplementary information of Janzen et al. 2018.
}
\usage{
time_error(J = NA, N = Inf, R = Inf, H_0 = 0.5, C = 1, t = 1, relative = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{J}{
      The number of junctions at time t
}
  \item{N}{
      Population Size
}
  \item{R}{
      Number of genetic markers
}
  \item{H_0}{
      Frequency of heterozygosity at t = 0
}
  \item{C}{
      Mean number of crossovers per meiosis (e.g. size in Morgan of the chromosome)
}
  \item{t}{
      Inferred time
}
  \item{relative}{
      Boolean flag, if TRUE: return the relative error, if FALSE: return error in generations
}
}

\value{
  the error in the time estimate
}
\author{
Thijs Janzen
}

\examples{
  time_error(J = 100, N = Inf, R = 1000, H_0 = 0.5, C = 1)
}

\keyword{analytic}
\keyword{time}
\keyword{time error}
