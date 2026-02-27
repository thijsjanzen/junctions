# Calculate the average number of junctions during backcrossing

Calculate the expected number of junctions after t generations, in a
backcrossing mating scheme.

## Usage

``` r
number_of_junctions_backcross(H_0 = 0.5, C = 1, t = 100)
```

## Arguments

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
jt <-  number_of_junctions_backcross(H_0 = 0.1, C = 1, t = 5)
```
