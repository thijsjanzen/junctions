# Estimate the error in the time estimate

Calculate the error in the estimate of the onset of hybridization,
following Equations 3 & 4 in the Supplementary information of Janzen et
al. 2018.

## Usage

``` r
time_error(t = NA, N = Inf, R = Inf, H_0 = 0.5, C = 1, relative = TRUE)
```

## Arguments

- t:

  Inferred time

- N:

  Population Size

- R:

  Number of genetic markers

- H_0:

  Frequency of heterozygosity at t = 0

- C:

  Mean number of crossovers per meiosis (e.g. size in Morgan of the
  chromosome)

- relative:

  Boolean flag, if TRUE: return the relative error, if FALSE: return
  error in generations

## Value

Expected error in the time estimate
