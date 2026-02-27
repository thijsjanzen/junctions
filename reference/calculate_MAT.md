# Function to calculate the maximum accurate time

Function that calculates the maximum time after hybridization after
which the number of junctions can still be reliably used to estimate the
onset of hybridization. This is following equation 15 in Janzen et al.
2018.

## Usage

``` r
calculate_mat(N = Inf, R = Inf, H_0 = 0.5, C = 1)
```

## Arguments

- N:

  Population Size

- R:

  Number of genetic markers

- H_0:

  Frequency of heterozygosity at t = 0

- C:

  Mean number of crossovers per meiosis (e.g. size in Morgan of the
  chromosome)

## Value

The maximum accurate time

## Examples

``` r
calculate_mat(N = Inf, R = 1000, H_0 = 0.5, C = 1)
#> [1] 6211.5
```
