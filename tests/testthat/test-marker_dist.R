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

  focal_time <- 500
  time_est <- c()

  for(r in seq_len(num_repl)) {
    sim_results <- sim_inf_chrom(pop_size = N,
                                 initial_heterozygosity = 0.5,
                                 total_runtime = total_runtime,
                                 morgan = 1,
                                 markers = R,
                                 seed = r)

    found_obs <- rbind(found_obs, sim_results$detectedJunctions)

    sim_markers <- sim_results$markers
    expected_junctions <-
      number_of_junctions_markers(N = N,
                                  H_0 = 0.5,
                                  C = 1,
                                  t = 0:total_runtime,
                                  marker_distribution = sim_markers)

    found_exp <- rbind(found_exp, expected_junctions)

    ers)
    time_est <- c(time_est, estimated_time)

    cat(r,"\n")
  }

  found_obs <- colMeans(found_obs)
  found_exp <- colMeans(found_exp)

  for(i in seq_along(found_obs)) {
    testthat::expect_equal(found_obs[i], found_exp[i], tolerance = 1)
  }

  testthat::expect_equal(mean(time_est), focal_time, tolerance = 1)
})

test_that("estimate time", {

  N <- 1e3
  total_runtime = 100


  sim_results <- sim_inf_chrom(pop_size = N,
                               initial_heterozygosity = 0.5,
                               total_runtime = total_runtime,
                               morgan = 1,
                               markers = R,
                               seed = 444)

  focal_j <- tail(sim_results$avgJunctions,1 )
  estimated_time <- estimate_time_markers(J = focal_j,
                                          N = N,
                                          H_0 = 0.5,
                                          C = 1,
                                          marker_distribution = sim_markers)

  testthat::expect_equal(estimated_time, total_runtime, tolerance = 1)
})