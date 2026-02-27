# Function to simulate data using a back crossing scheme

Individual based simulation of the accumulation of junctions, under a
back crossing scheme

## Usage

``` r
sim_backcrossing(
  population_size = 100,
  freq_ancestor_1 = 0.5,
  total_runtime = 5,
  size_in_morgan = 1,
  number_of_markers = 100,
  seed = 6,
  time_points = -1
)
```

## Arguments

- population_size:

  Population size

- freq_ancestor_1:

  Frequency of ancestor 1 at t = 0

- total_runtime:

  Number of generations to simulate

- size_in_morgan:

  Mean number of crossovers per meiosis (e.g. size in Morgan of the
  chromosome)

- number_of_markers:

  number of molecular markers

- seed:

  Seed of the pseudo-random number generator

- time_points:

  vector with time points at which local ancestry has to be recorded to
  be returned at the end of the simulation. If left at -1, ancestry is
  recorded at every generation (computationally heavy).

## Value

List with five entries: average_junctions: average number of junctions
over time, detected_junctions: average number of detected junctions,
given the markers. markers: vector with the locations of the molecular
markers, junction_distribution: distribution of junctions per time step
average_heterozygosity: average heterozygosity.

## Examples

``` r
sim_backcrossing(population_size = 100,
                       total_runtime = 5,
                       size_in_morgan = 1,
                       number_of_markers = 100,
                       seed = 6,
                       time_points = 1:5)
#> $average_junctions
#> [1] 0.485 0.440 0.310 0.180
#> 
#> $detected_junctions
#> [1] 0.465 0.415 0.295 0.170
#> 
#> $markers
#>   [1] 0.01180642 0.01933798 0.03497799 0.04376882 0.04464693 0.06775731
#>   [7] 0.08209470 0.10351912 0.11354000 0.12224747 0.12887618 0.13357049
#>  [13] 0.13790556 0.14434618 0.14638019 0.15897425 0.16708516 0.16776536
#>  [19] 0.17159321 0.19303798 0.19889418 0.20688505 0.20875709 0.22089998
#>  [25] 0.23957689 0.25365198 0.25883235 0.26029875 0.27446782 0.27980986
#>  [31] 0.34747154 0.34951547 0.35204047 0.37812669 0.38477364 0.39630588
#>  [37] 0.42294457 0.42592437 0.43387928 0.44607088 0.44837856 0.44922129
#>  [43] 0.45086201 0.47754811 0.48334898 0.49151114 0.50652364 0.51913075
#>  [49] 0.52846650 0.53509150 0.54148667 0.54480760 0.55031412 0.56009784
#>  [55] 0.56861911 0.57127711 0.58964303 0.59084356 0.59543129 0.63511870
#>  [61] 0.63747898 0.63927276 0.64007345 0.64138804 0.64373512 0.64495380
#>  [67] 0.65056443 0.65096715 0.65745963 0.66013086 0.67010347 0.67782402
#>  [73] 0.68643371 0.68887698 0.70851656 0.73601246 0.73717091 0.73777945
#>  [79] 0.77606555 0.78322536 0.78563072 0.80229424 0.80282633 0.80720381
#>  [85] 0.81155567 0.82451640 0.84499916 0.85470518 0.85704208 0.86058915
#>  [91] 0.86325870 0.87209400 0.88216151 0.91082443 0.91456731 0.91797463
#>  [97] 0.93506913 0.94177503 0.98056468 0.99789406
#> 
#> $junction_distribution
#> $junction_distribution[[1]]
#>   [1] 1 0 2 1 1 4 0 1 3 4 0 0 0 1 1 2 2 2 0 1 1 0 2 1 1 1 0 0 1 0 0 2 2 1 1 2 1
#>  [38] 0 2 1 0 2 0 0 1 0 0 0 1 0 1 0 1 0 1 1 2 0 2 0 3 1 1 1 0 0 1 0 3 1 2 2 0 1
#>  [75] 0 0 2 2 0 3 1 1 2 0 0 0 0 0 4 0 1 0 0 0 2 0 1 0 0 1
#> 
#> $junction_distribution[[2]]
#>   [1] 0 0 2 3 3 0 0 1 0 3 0 0 0 0 2 0 2 2 2 0 2 1 0 1 1 1 1 0 0 1 2 3 0 2 0 4 0
#>  [38] 0 0 1 1 2 1 0 1 1 2 2 1 2 0 1 0 2 1 2 0 1 0 0 1 2 1 0 0 1 0 0 0 1 1 0 2 0
#>  [75] 0 0 2 1 1 0 2 1 1 0 2 1 2 0 0 0 0 1 0 0 0 0 0 0 0 0
#> 
#> $junction_distribution[[3]]
#>   [1] 0 0 1 0 2 4 0 3 0 0 0 0 0 1 0 0 2 0 0 2 0 0 0 0 0 0 2 0 0 1 0 1 0 1 0 0 0
#>  [38] 0 0 0 2 0 1 0 0 0 0 2 0 0 0 0 0 0 1 0 0 0 0 1 1 0 2 2 0 0 2 1 0 0 2 3 1 0
#>  [75] 2 1 0 0 0 0 2 0 0 1 2 0 0 2 1 3 0 1 0 0 0 1 2 0 0 0
#> 
#> $junction_distribution[[4]]
#>   [1] 0 0 0 1 0 0 0 0 1 0 0 1 1 0 0 0 2 0 0 0 0 2 0 2 1 0 1 0 0 2 0 0 0 0 0 0 0
#>  [38] 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 1 0 0 0 2 0 0 0 0
#>  [75] 0 1 2 0 0 0 1 0 0 0 2 0 1 0 0 0 0 2 0 2 2 0 0 0 0 0
#> 
#> 
#> $average_heterozygosity
#> [1] 0.4698 0.2217 0.1111 0.0450
#> 
```
