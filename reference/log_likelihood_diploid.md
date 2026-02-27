# calculate the log likelihood of observing diploid ancestry data.

Calculates the log likelihood of observing the phased data, given the
population size, initial heterozygosity and time since admixture

## Usage

``` r
log_likelihood_diploid(
  local_anc_matrix,
  pop_size,
  freq_ancestor_1 = 0.5,
  t,
  phased = FALSE,
  num_threads = 1
)
```

## Arguments

- local_anc_matrix:

  a matrix with four columns: column 1) chromosome indicator, 2)
  location of marker in Morgan on respective chromosome 3) ancestry at
  chromosome 4) ancestry at chromosome 2.

- pop_size:

  population size

- freq_ancestor_1:

  Frequency of ancestor 1 at t = 0

- t:

  time since admixture

- phased:

  is the data phased or not? default is false.

- num_threads:

  number of threads, default is one thread. Set to -1 to use all
  available threads.

## Value

log likelihood
