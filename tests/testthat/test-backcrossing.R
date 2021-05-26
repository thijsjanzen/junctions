context("backcross")
test_that("backcross, use", {
  testthat::skip_on_os("solaris")
  vx <- sim_backcrossing(population_size = 10000,
                           total_runtime = 10,
                           freq_ancestor_1 = 0.5,
                           time_points = 1:10,
                           seed = 42)
  t <- 1:9
  expected_heterozygosity <-  2 ^ (-t)

  testthat::expect_equal(expected_heterozygosity,
                         vx$average_heterozygosity, tolerance = 0.01)

  expected_junctions <- number_of_junctions_backcross(t = 1:9)
  observed_junctions <- vx$average_junctions

  testthat::expect_equal(expected_junctions,
                         observed_junctions, tolerance = 0.05)


  vx <- sim_backcrossing(population_size = 10000,
                         total_runtime = 10,
                         freq_ancestor_1 = 0.2,
                         seed = 42)
})
