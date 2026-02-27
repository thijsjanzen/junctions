# Calculate the limit of the number of junctions

Calculate the average number of junctions after an infinite number of
generations, provided information on the initial heterozygosity,
population size and the number of generations.

## Usage

``` r
calc_k(N = Inf, R = Inf, H_0 = 0.5, C = 1)
```

## Arguments

- N:

  population size

- R:

  number of markers

- H_0:

  initial heterozygosity (at the time of admixture)

- C:

  Mean number of crossovers per meiosis (e.g. size in Morgan of the
  chromosome)

## Value

The number of junctions for at time = infinity

## Examples

``` r
k <-  calc_k(N = 100, R = 1000, H_0 = 0.5, C = 1)
```
