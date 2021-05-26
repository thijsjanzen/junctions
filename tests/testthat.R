library(testthat)
library(junctions)

backend <- Sys.getenv("RCPP_PARALLEL_BACKEND", unset = NA)
if (is.na(backend))
  Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

writeLines(paste("Using backend:", Sys.getenv("RCPP_PARALLEL_BACKEND")))

test_check("junctions")
