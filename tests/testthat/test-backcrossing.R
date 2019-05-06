context("backcross")
test_that("backcross, use", {

  vx <- sim_backcrossing(population_size = 10000,
                              max_time = 10,
                              freq_ancestor_1 = 0.5,
                              time_points = 1:10)

  t <- 1:10
  expected_heterozygosity <- 0.5 * t*  2^(-t)
  testthat::expect_equal(expected_heterozygosity,
                         vx$average_heterozygosity)
})
