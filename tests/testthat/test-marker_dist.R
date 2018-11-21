context("marker_distribution")
test_that("calculate_J and error", {

  R <- 1000
  e_j <- number_of_junctions(R = R,
                             t = 200)

  markers <- seq(0,1, length.out = R+1)
  d_e_j <- number_of_junctions_markers(t = 200,
                                       marker_distribution = markers)

  testthat::expect_equal(e_j, d_e_j)

  num_repl <- 100

  N <- 50
  R <- 1000
  total_runtime = 1000

  found_obs <- c()
  found_exp <- c()
  for(r in seq_len(num_repl)) {
    sim_results <- sim_inf_chrom(pop_size = N,
                                 init_heterozygosity = 0.5,
                                 run_time = total_runtime,
                                 size_in_Morgan = 1,
                                 markers = R)

    found_obs <- rbind(found_obs, sim_results$detectedJunctions)

    sim_markers <- sim_results$markers
    expected_junctions <-
      number_of_junctions_markers(N = N,
                                  H_0 = 0.5,
                                  C = 1,
                                  t = total_runtime,
                                  marker_distribution = sim_markers)

    found_exp <- rbind(found_exp, expected_junctions)
  }

  found_obs <- colMeans(found_obs)
  found_exp <- colMeans(found_exp)

  for(i in seq_along(found_obs)) {
    testthat::expect_equal(found_obs[i], found_exp[i])
  }
})