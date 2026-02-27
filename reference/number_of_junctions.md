# Calculate the average number of junctions

Calculate the average number of junctions in a single chromosome after t
generations, provided information on the initial heterozygosity,
population size and the number of generations.

## Usage

``` r
number_of_junctions(N = Inf, R = Inf, H_0 = 0.5, C = 1, t = 100)
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

- t:

  Time since admixture

## Value

Estimated number of junctions at time t

## Examples

``` r
jt <-  number_of_junctions(N = 100, R = 1000, H_0 = 0.5, C = 1, t = 1000)
jt2 <- number_of_junctions(N = 100, R = 1000, H_0 = 0.5, C = 1, t = 0:1000)
```
