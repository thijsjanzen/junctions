context("marker_distribution")
test_that("markers, use", {
  testthat::skip_on_os("solaris")

  num_markers <- 1000
  e_j <- number_of_junctions(R = num_markers,
                             t = 200)

  markers <- seq(0, 1, length.out = num_markers + 1)
  d_e_j <- number_of_junctions_markers(t = 200,
                                       marker_distribution = markers)

  testthat::expect_equal(e_j, d_e_j)

  num_repl <- 10

  N <- 1000   # nolint
  num_markers <- 1000
  total_runtime <- 20

  found_obs <- c()
  found_exp <- c()

  for (r in seq_len(num_repl)) {
    sim_results <- sim_inf_chrom(pop_size = N,
                                 freq_ancestor_1 = 0.5,
                                 total_runtime = total_runtime,
                                 morgan = 1,
                                 markers = num_markers,
                                 seed = r)

    found_obs <- rbind(found_obs, sim_results$detectedJunctions) # nolint

    sim_markers <- sim_results$markers
    expected_junctions <-
      number_of_junctions_markers(N = N,
                                  H_0 = 0.5,
                                  t = 0:total_runtime,
                                  marker_distribution = sim_markers)

    found_exp <- rbind(found_exp, expected_junctions)
  }

  found_obs <- colMeans(found_obs)
  found_exp <- colMeans(found_exp)

  for (i in seq_along(found_obs)) {
    testthat::expect_equal(found_obs[i], found_exp[i], tolerance = 0.1)
  }
})

test_that("markers, abuse", {
  testthat::expect_error(
    number_of_junctions_markers(marker_distribution = 0.5)
  )
})

test_that("marker dist", {
  testthat::skip_on_os("solaris")
  message("marker dist")
  vx <- sim_phased_unphased(total_runtime = 2,
                            markers = 10,
                            num_indiv_sampled = 1,
                            time_points = 1)

  testthat::expect_equal(length(vx$location), 10)
  diff_loc <- diff(vx$location)
  testthat::expect_false(diff_loc[1] == diff_loc[2])

  vy <- sim_phased_unphased(total_runtime = 2,
                            markers = -9,
                            num_indiv_sampled = 1,
                            time_points = 1)

  testthat::expect_equal(length(vy$location), 9)
  diff_loc <- diff(vy$location)
  testthat::expect_true(diff_loc[1] == diff_loc[2])

  a <- get_num_markers(-9)
  b <- diff(a)
  testthat::expect_true(all.equal(diff_loc, b))
})
