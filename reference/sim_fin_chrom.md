# Individual Based Simulation of the accumulation of junctions

Individual based simulation of the accumulation of junctions for a
chromosome with regularly distributed markers.

## Usage

``` r
sim_fin_chrom(
  pop_size = 100,
  freq_ancestor_1 = 0.5,
  total_runtime = 100,
  morgan = 1,
  seed = 42,
  R = 100
)
```

## Arguments

- pop_size:

  Population Size

- freq_ancestor_1:

  Frequency of ancestor 1 at t = 0

- total_runtime:

  Maximum time after which the simulation is to be stopped

- morgan:

  Mean number of crossovers per meiosis (e.g. size in Morgan of the
  chromosome)

- seed:

  Seed of the pseudo-random number generator

- R:

  Number of regularly distributed markers

## Value

- avgJunctions:

  vector of the average number of junctions at time = \[0,
  total_runtime\]

## Examples

``` r
sim_fin_chrom(pop_size = 100, freq_ancestor_1 = 0.5,
                   total_runtime = 10, morgan = 1, seed = 42,
                   R = 100)
#> $avgJunctions
#>  [1] 0.000 0.405 0.925 1.240 1.655 2.140 2.590 3.325 3.915 4.365 4.660
#> 
```
