# junctions <img src="pics/junctions_sticker3.png" align="right" width="180" />
Individual based simulations of hybridizing populations, where the accumulation of junctions is tracked. Furthermore, mathematical equations are provided to verify simulation outcomes. Both simulations and mathematical equations are based on Janzen et al. (2018) 

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/junctions)](https://cran.r-project.org/package=junctions)
[![](http://cranlogs.r-pkg.org/badges/grand-total/junctions)](https://cran.r-project.org/package=junctions)
[![](http://cranlogs.r-pkg.org/badges/junctions)](https://cran.r-project.org/package=junctions)

Branch|[![Travis CI logo](pics/TravisCI.png)](https://travis-ci.org)|[![AppVeyor logo](pics/AppVeyor.png)](https://www.appveyor.com)|[![Codecov logo](pics/Codecov.png)](https://www.codecov.io)
---|---|---|---
master|[![Build Status](https://travis-ci.org/thijsjanzen/junctions.svg?branch=master)](https://travis-ci.org/thijsjanzen/junctions)|[![Build status](https://ci.appveyor.com/api/projects/status/rt9856tv3pi87sms?svg=true)](https://ci.appveyor.com/project/thijsjanzen/junctions)|[![codecov.io](https://codecov.io/gh/thijsjanzen/junctions/branch/master/graph/badge.svg)](https://codecov.io/gh/thijsjanzen/junctions)




## references
Janzen, T. , Nolte, A. W. and Traulsen, A. (2018), The breakdown of genomic ancestry blocks in hybrid lineages given a finite number of recombination sites. Evolution, 72: 735-750. https://doi.org/10.1111/evo.13436

Lavretsky, P, Janzen, T. and McCracken, KG.  (2019) Identifying hybrids & the genomics of hybridization: Mallards & American black ducks of Eastern North America. Ecology and Evolution 9: 3470-3490. https://doi.org/10.1002/ece3.4981

## Updates
Version 1.6 :  Improved the recombination function run twice as fast <br />
Version 1.5.1: Added option to track the true number of junctions <br />
Version 1.5  : Added simulation functions to simulate phased an unphased data, including phasing error <br />
Version 1.5  : Added support for inferring the time since admixture based on phased and unphased data. <br />
Version 1.4  : Added support for estimating the number of junctions, and simulating the number of junctions, under a backcrossing scheme, using the code supplied in Lavretsky et al. 2019. <br />
Version 1.3  : Added support for estimating the time since admixture using unphased data. <br />
Version 1.3  : Added individual based simulations returning phased and unphased data. <br />
Version 1.3  : Updated entire package to Roxygen. <br />
Version 1.2  : Added support for estimating the expected number of junctions for arbitrarily distributed markers. <br />
Version 1.1  : Updated random number generation for picking recombination sites. Previous implementation was limited to 6 digit precision, current precision is at least double that, minimizing the probability of recombination occur twice in the same location for an infinite chromosome. <br />
