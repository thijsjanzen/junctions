# log likelihood of the time since admixture for a haploid genome

log likelihood of the time since admixture for a set of single
chromosomes (for ex. in Yeast).

## Usage

``` r
log_likelihood_haploid(ancestry_matrix, N = 1000, freq_ancestor_1 = 0.5, t = 2)
```

## Arguments

- ancestry_matrix:

  matrix with 3 columns, column 1 = chromosome, column 2 = location in
  Morgan, column 3 = ancestry.

- N:

  Population Size

- freq_ancestor_1:

  Frequency of ancestor 1 at t = 0

- t:

  time since admixture

## Value

loglikelihood
