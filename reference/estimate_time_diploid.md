# estimates the time since admixture, given diploid ancestry data.

Calculates the time since admixture, given unphased ancestry data.

## Usage

``` r
estimate_time_diploid(
  ancestry_information,
  analysis_type = "individuals",
  phased = FALSE,
  pop_size = 1000,
  freq_ancestor_1 = 0.5,
  lower_lim = 2,
  upper_lim = 2000,
  num_threads = 1,
  verbose = FALSE
)
```

## Arguments

- ancestry_information:

  a matrix with five columns: column 1) indicator of individual,
  column 2) indicator of chromosome, 3) location of marker in Morgan, 4)
  ancestry at chromosome 5) ancestry at chromosome 2.

- analysis_type:

  how should the data be broken down? there are multiple options:
  "individuals" - time is inferred for each individual separately,
  grouping all chromosomes together that belong to the same individual.
  "chromosomes" - time is inferred for each chromosome separately,
  grouping chromosomes together belonging from separate individuals.
  "separate" - time is inferred for each chromosome from each individual
  separately, "all" - time is inferred jointly for all chromosomes and
  individuals, grouping all chromosomes and individuals together.

- phased:

  is the data phased?

- pop_size:

  population size

- freq_ancestor_1:

  Frequency of ancestor 1 at t = 0

- lower_lim:

  lower limit of the optimization algorithm. Increase if the expected
  admixture time is relatively ancient

- upper_lim:

  upper limit of hte optimization algorithm. If set too large, recent
  admixture events can be overlooked - best to set as low as possible.

- num_threads:

  num_threads, default is all threads. 5 threads is recommended.

- verbose:

  display intermediate output? Default = FALSE
