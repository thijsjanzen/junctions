context("unphased")
test_that("unphased, use", {

  vx <- sim_phased_unphased(pop_size = 100,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 201,
                            size_in_morgan = 1,
                            number_of_markers = 1000,
                            time_points = c(100, 200),
                            seed = 42)

  num_indiv <- length(unique(vx$results[,2]))
  testthat::expect_equal(num_indiv, 10)
  testthat::expect_equal(length(unique(vx$results[,1])), 2)

  found <- c()
  focal_time <- 100
  for(indiv in 0:9) {
    local_data <- subset(vx$results, vx$results[,2] == indiv &
                           vx$results[,1] == focal_time)
    estim_time <- estimate_time_unphased(local_data[,4:5], local_data[,3],
                                         pop_size = 100, freq_ancestor_1 = 0.5,
                                         max_t = 200)
    found[indiv+1] <- estim_time$minimum
  }

  testthat::expect_equal(mean(found), focal_time, tolerance = 10)

  local_data <- subset(vx$results, vx$results[,2] == 0 &
                         vx$results[,1] == 100)
  ll_100 <- unphased_log_likelihood(local_data[,4:5], local_data[,3],
                                    pop_size = 100,
                                    freq_ancestor_1 = 0.5,
                                    t = 100)

  ll_200 <- unphased_log_likelihood(local_data[,4:5], local_data[,3],
                                    pop_size = 100,
                                    freq_ancestor_1 = 0.5,
                                    t = 200)

  testthat::expect_gte(ll_100, ll_200)
})