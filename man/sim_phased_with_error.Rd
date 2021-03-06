% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_phased_with_error.R
\name{sim_phased_with_error}
\alias{sim_phased_with_error}
\title{Individual Based Simulation of the accumulation of junctions, given phasing
error}
\usage{
sim_phased_with_error(
  pop_size = 100,
  freq_ancestor_1 = 0.5,
  total_runtime = 100,
  size_in_morgan = 1,
  markers = 100,
  time_points = -1,
  seed = NULL,
  num_threads = 1,
  verbose = FALSE,
  num_indiv_sampled = 10,
  coverage = 1,
  error_rate = 0
)
}
\arguments{
\item{pop_size}{Population Size}

\item{freq_ancestor_1}{Frequency of ancestor 1 at t = 0}

\item{total_runtime}{Maximum time after which the simulation is to be stopped}

\item{size_in_morgan}{Mean number of crossovers per meiosis (e.g. size in
Morgan of the chromosome)}

\item{markers}{If a single number is provided, the number is used as the
total number of markers generated either randomly, or using a regular
distribution (a regular distribution is chosen if the number is negative). If
a vector is provided, that vector is used.}

\item{time_points}{vector with time points at which local ancestry has to be
recorded to be returned at the end of the simulation. If left at -1,
ancestry is recorded at every generation (computationally heavy).}

\item{seed}{Seed of the pseudo-random number generator}

\item{num_threads}{number of threads, deafault is 1. -1 takes all available
threads.}

\item{verbose}{displays a progress bar, default is FALSE}

\item{num_indiv_sampled}{the number of individuals sampled at each time point
to be genotyped}

\item{coverage}{fraction of markers that can be succesfully phased}

\item{error_rate}{fraction of markers that are erroneously
phased (e.g. swapped)}
}
\value{
a tibble with five columns: [time, individual, marker location,
                                ancestry chromosome 1, ancestry chromosome 2]
}
\description{
Individual based simulation of the accumulation of junctions,
returning phased and unphased data. Phased data is generated using a
user-provided error rate. Ancestry on both chromosomes of 10 randomly sampled
individuals per generations is returned.
}
\examples{
\dontrun{
sim_phased_with_error(pop_size = 100, freq_ancestor_1 = 0.5,
                    total_runtime = 1000, size_in_morgan = 1,
                    markers = 10, time_points = c(100),
                    seed = 42, coverage = 0.9, error_rate = 0.01)
 }
}
