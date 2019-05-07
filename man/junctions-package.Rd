\name{junctions-package}
\alias{junctions-package}
\alias{junctions}
\docType{package}
\title{
  Extending The Theory of Junctions
}
\description{
  The theory of junctions is extended by this package by including the effect of a finite number of recombination sites along the chromosome. The package provides functions to calculate the estimated number of junctions, depending on the time since the onset of hybridization, population size, number of recombination sites, initial heterozygosity and the number of crossovers per meiosis.
}
\details{
  This package provides individual based simulations in order to simulate the accumulation of junctions over time, both for chromosomes with a finite and an infinite number of recombination sites. Furthermore, the package provides mathetematical tools to verify the outcomes of the individual based simulations.

  Update version 1.4: added support for estimating the number of junctions, and simulating the number of junctions, under a backcrossing scheme, using the code supplied in Lavretsky et al. 2019. \cr

  Update version 1.3: added support for estimating the time since admixture using unphased data. \cr
  Update version 1.3: added individual based simulations returning phased and unphased data. \cr
  Update version 1.3: Updated entire package to Roxygen. \cr

  Update version 1.2: added support for estimating the expected number of junctions for arbitrarily distributed markers. \cr

  Update version 1.1: updated underlying random number generator for picking recombination sites. The previous generator had limited precision, which could generate duplicate recombination sites. This update fixes that \cr
}
\author{
Maintainer: Thijs Janzen <thijsjanzen@gmail.com>
}
\references{
  Janzen, T. , Nolte, A. W. and Traulsen, A. (2018), The breakdown of genomic ancestry blocks in hybrid lineages given a finite number of recombination sites. Evolution, 72: 735-750. doi:10.1111/evo.13436 \cr

  Lavretsky, P, Janzen, T. and McCracken, KG.  (2019) Identifying hybrids & the genomics of hybridization: Mallards & American black ducks of Eastern North America. Ecology and Evolution 9: 3470-3490. doi:10.1002/ece3.4981 \cr

}
\keyword{ junctions }

