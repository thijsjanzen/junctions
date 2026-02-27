# estimate time using likelihood for a single chromosome

Estimate the time since the onset of hybridization, for a haploid genome

## Usage

``` r
estimate_time_haploid(
  ancestry_matrix,
  N = 1000,
  freq_ancestor_1 = 0.5,
  lower_lim = 2,
  upper_lim = 1000,
  verbose = FALSE
)
```

## Arguments

- ancestry_matrix:

  matrix with 3 columns, column 1 = chromosome, column 2 = location in
  Morgan, column 3 = ancestry.

- N:

  Population Size

- freq_ancestor_1:

  Frequency of ancestor 1 at t = 0

- lower_lim:

  lower limit of the optimization algorithm. Increase if the expected
  admixture time is relatively ancient

- upper_lim:

  upper limit of the optimization algorithm. If set too large, recent
  admixture events can be overlooked - best to set as low as possible.

- verbose:

  return verbose output

## Value

The number of generations passed since the onset of hybridization
