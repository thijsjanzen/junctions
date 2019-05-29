context("unphased")
test_that("unphased, use", {

  vx <- sim_phased_unphased(pop_size = 100,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 201,
                            size_in_morgan = 1,
                            number_of_markers = 1000,
                            time_points = c(100, 200),
                            seed = 42)

  num_indiv <- length(unique(vx$individual))
  testthat::expect_equal(num_indiv, 10)
  testthat::expect_equal(length(unique(vx$time)), 2)

  found <- c()
  focal_time <- 100
  for(indiv in 0:9) {
    local_data <- subset(vx, vx$individual == indiv &
                           vx$time == focal_time)


    estim_time <- estimate_time_unphased(cbind(local_data$anc_chrom_1, local_data$anc_chrom_2),
                                         local_data$location,
                                         pop_size = 100, freq_ancestor_1 = 0.5,
                                         upper_lim = 200)
    found[indiv+1] <- estim_time$minimum
  }

  testthat::expect_equal(mean(found), focal_time, tolerance = 10)



  local_data <- subset(vx, vx$individual == 0 &
                         vx$time == 100)
  ll_100 <- unphased_log_likelihood(cbind(local_data$anc_chrom_1, local_data$anc_chrom_2),
                                    local_data$location,
                                    pop_size = 100,
                                    freq_ancestor_1 = 0.5,
                                    t = 100)

  ll_200 <- unphased_log_likelihood(cbind(local_data$anc_chrom_1, local_data$anc_chrom_2),
                                    local_data$location,
                                    pop_size = 100,
                                    freq_ancestor_1 = 0.5,
                                    t = 200)

  testthat::expect_gte(ll_100, ll_200)

  vx <- sim_phased_unphased(pop_size = 1000,
                            freq_ancestor_1 = 0.1,
                            total_runtime = 30,
                            size_in_morgan = 1,
                            number_of_markers = 1000,
                            time_points = c(30),
                            seed = 421)

  local_data <- subset(vx, vx$individual == 0 &
                         vx$time == 30)

  ll_30 <- unphased_log_likelihood(cbind(local_data$anc_chrom_1, local_data$anc_chrom_2),
                                    local_data$location,
                                    pop_size = 1000,
                                    freq_ancestor_1 = 0.1,
                                    t = 30)

  ll_60 <- unphased_log_likelihood(cbind(local_data$anc_chrom_1, local_data$anc_chrom_2),
                                    local_data$location,
                                    pop_size = 1000,
                                    freq_ancestor_1 = 0.1,
                                    t = 60)
  testthat::expect_gte(ll_30, ll_60)

  ll_inf <- unphased_log_likelihood(cbind(local_data$anc_chrom_1, local_data$anc_chrom_2),
                                    local_data$location,
                                    pop_size = 1000,
                                    freq_ancestor_1 = 0.1,
                                    t = 0)
  testthat::expect_true(is.infinite(ll_inf))
})