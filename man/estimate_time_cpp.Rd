% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{estimate_time_cpp}
\alias{estimate_time_cpp}
\title{function to calculate log likelihood using cpp}
\usage{
estimate_time_cpp(
  local_anc_matrix,
  locations,
  pop_size,
  freq_ancestor_1,
  lower_lim,
  upper_lim,
  verbose,
  phased,
  num_threads = -1L
)
}
\arguments{
\item{local_anc_matrix}{local ancestry matrix}

\item{locations}{locations of markers}

\item{pop_size}{population size}

\item{freq_ancestor_1}{frequency of the most common ancestor}

\item{lower_lim}{lower limit}

\item{upper_lim}{upper limit}

\item{verbose}{use verbose output}

\item{phased}{is the data phased or not?}

\item{num_threads, }{default is all threads. 5 threads is recommended.}
}
\description{
function to calculate log likelihood using cpp
}
