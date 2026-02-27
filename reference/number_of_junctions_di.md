# Calculate the expected number of junctions between two markers separated by a given amount of recombination

Calculate the expected number of junctions after t generations, provided
information on the initial heterozygosity, population size, the number
of generations since the onset of admixture and the distance between two
markers.

## Usage

``` r
number_of_junctions_di(N = Inf, H_0 = 0.5, t = 100, di = 1e-06)
```

## Arguments

- N:

  Population Size

- H_0:

  Frequency of heterozygosity at t = 0

- t:

  Time since admixture

- di:

  Distance between two markers in Morgan

## Value

Estimated number of junctions at time t

## Examples

``` r
number_of_junctions_di(N = 100, H_0 = 0.5, t = 1000, di = 0.01)
#> [1] 0.3333332
```
