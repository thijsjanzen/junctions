# Estimate the time since the onset of hybridization, using the number of junctions

Estimate the time since the onset of hybridization, following equation
14 in Janzen et al. 2018

## Usage

``` r
estimate_time(J = NA, N = Inf, R = Inf, H_0 = 0.5, C = 1)
```

## Arguments

- J:

  The observed number of junctions

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

The number of generations passed since the onset of hybridization

## Examples

``` r
J <- number_of_junctions(N = 100, R = 1000, H_0 = 0.5, C = 1, t = 200)
estimate_time(J = J, N = 100, R = 1000, H_0 = 0.5, C = 1)
#> [1] 200
# should be 200 again
```
