This package provides individual based simulations in order to simulate the
accumulation of junctions over time, both for chromosomes with a finite and an
infinite number of recombination sites. Furthermore, the package provides
mathematical tools to verify the outcomes of the individual based simulations.

### version 2.1.0
- updated tbb::task_scheduler_init to tbb::global_control 
- improved linting of cpp code
- sped up tests
- added CITATION file
- added NEWS file

### version 2.0.2
- simplified some tests  

### version 2.0
- merged many functions with similar functionality
- added vignette that provides overview of all functionality.  

### version 1.9
- added c++ versions of the unphased and phased likelihoods.  

### version 1.8
- added multithreading using the TBB / RcppParallel library.  

### version 1.7
- further improved the recombination function following Hanno Hildenbrandt's 
suggestions  

### version 1.6
- improved the recombination function to run twice as fast  

### version 1.5.1
- added option to track the true number of junctions  

### version 1.5
- added support for inferring the time since admixture based on phased and
unphased data. 
- included simulation functions to simulate appropriate data (e.g. phased and
unphased).  

### version 1.4
- added support for estimating the number of junctions under a backcrossing
scheme
- added support for simulating the number of junctions under a backcrossing
scheme
- both based on the code supplied in Lavretsky et al. 2019.  

### version 1.3
- added support for estimating the time since admixture using unphased data.  
- added individual based simulations returning phased and unphased data.  
- Updated entire package to Roxygen.  

### version 1.2
- added support for estimating the expected number of junctions for
arbitrarily distributed markers.  

### version 1.1 
- updated underlying random number generator for picking recombination sites. 
The previous generator had limited precision, which could generate duplicate
recombination sites. This update fixes that.