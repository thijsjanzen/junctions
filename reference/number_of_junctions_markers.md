# Calculate the expected total number of junctions in a chromosome, given the distribution of markers

Calculate the expected number of junctions after t generations, provided
information on the initial heterozygosity, population size, the number
of generations since the onset of admixture and the distribution of
markers.

## Usage

``` r
number_of_junctions_markers(
  N = Inf,
  H_0 = 0.5,
  t = 100,
  marker_distribution = NA
)
```

## Arguments

- N:

  Population Size

- H_0:

  Frequency of heterozygosity at t = 0

- t:

  Time since admixture

- marker_distribution:

  A vector containing the position of all markers in Morgan.

## Value

Estimated number of observed junctions at time t

## Examples

``` r
markers <- seq(from = 0, to = 1, length.out = 1000)
jt <-  number_of_junctions_markers(N = 100,
                                  H_0 = 0.5,
                                  t = 1000,
                                  marker_distribution = markers)
random_markers <- sort(runif(1000, 0, 1))
jt2 <- number_of_junctions_markers(N = 100,
                                  H_0 = 0.5,
                                  t = 1000,
                                  marker_distribution = random_markers)
```
