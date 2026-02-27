# Individual Based Simulation of the accumulation of junctions

Individual based simulation of the accumulation of junctions, returning
phased and unphased data. Ancestry on both chromosomes of 10 randomly
sampled individuals per generations is returned.

## Usage

``` r
sim_phased_unphased(
  pop_size = 100,
  freq_ancestor_1 = 0.5,
  total_runtime = 100,
  size_in_morgan = 1,
  markers = 100,
  time_points = -1,
  num_threads = 1,
  verbose = FALSE,
  record_true_junctions = FALSE,
  num_indiv_sampled = 10,
  coverage = 1,
  error_rate = 0
)
```

## Arguments

- pop_size:

  Population Size

- freq_ancestor_1:

  Frequency of ancestor 1 at t = 0

- total_runtime:

  Maximum time after which the simulation is to be stopped

- size_in_morgan:

  Mean number of crossovers per meiosis (e.g. size in Morgan of the
  chromosome)

- markers:

  If a single number is provided, the number is used as the total number
  of markers generated either randomly, or using a regular distribution
  (a regular distribution is chosen if the number is negative). If a
  vector is provided, that vector is used.

- time_points:

  vector with time points at which local ancestry has to be recorded to
  be returned at the end of the simulation. If left at -1, ancestry is
  recorded at every generation (computationally heavy).

- num_threads:

  default is 1. -1 takes all available threads.

- verbose:

  displays a progress bar

- record_true_junctions:

  keep track of the true number of junctions?

- num_indiv_sampled:

  the number of individuals sampled at each time point to be genotyped

- coverage:

  fraction of markers that can be succesfully phased

- error_rate:

  fraction of markers that are erroneously phased (e.g. swapped)

## Value

a tibble with five columns: \[time, individual, marker location,
ancestry chromosome 1, ancestry chromosome 2\]

## Examples

``` r
if (FALSE) { # \dontrun{
sim_phased_unphased(pop_size = 100, freq_ancestor_1 = 0.5,
                    total_runtime = 10, size_in_morgan = 1,
                    markers = 10, time_points = c(0, 5, 10),
                    num_threads = 1)
} # }
```
