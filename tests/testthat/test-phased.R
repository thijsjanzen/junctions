context("phased")
test_that("phased, use", {
  testthat::skip_on_os("solaris")
  vx <- sim_phased_unphased(pop_size = 100,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 120,
                            size_in_morgan = 1,
                            markers = 1000,
                            time_points = c(100, 120))

  num_indiv <- length(unique(vx$individual))
  testthat::expect_equal(num_indiv, 10)
  testthat::expect_equal(length(unique(vx$time)), 2)

  found <- c()
  focal_time <- 100
  for (indiv in 0:9) {
    local_data <- subset(vx, vx$individual == indiv &
                           vx$time == focal_time)


    estim_time <- estimate_time_diploid(cbind(1, 1,
                                              local_data$location,
                                              local_data$anc_chrom_1,
                                             local_data$anc_chrom_2),
                                       phased = TRUE,
                                       pop_size = 100,
                                       freq_ancestor_1 = 0.5,
                                       upper_lim = 200)
    found[indiv + 1] <- estim_time$time[[1]]
  }

  testthat::expect_equal(mean(found), focal_time, tolerance = 10)


  local_data <- subset(vx, vx$individual == 0 &
                         vx$time == 100)
  ll_100 <- log_likelihood_diploid(cbind(1,
                                         local_data$location,
                                         local_data$anc_chrom_1,
                                       local_data$anc_chrom_2),
                                 phased = TRUE,
                                 pop_size = 100,
                                 freq_ancestor_1 = 0.5,
                                 t = 100)

  ll_2000 <- log_likelihood_diploid(cbind(1,
                                         local_data$location,
                                         local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                  phased = TRUE,
                                  pop_size = 100,
                                  freq_ancestor_1 = 0.5,
                                 t = 2000)

  testthat::expect_gte(ll_100, ll_2000)

  vx <- sim_phased_unphased(pop_size = 1000,
                            freq_ancestor_1 = 0.1,
                            total_runtime = 30,
                            size_in_morgan = 1,
                            markers = 1000,
                            time_points = c(30))

  local_data <- subset(vx, vx$individual == 0 &
                         vx$time == 30)

  ll_30 <- log_likelihood_diploid(cbind(1,
                                        local_data$location,
                                        local_data$anc_chrom_1,
                                        local_data$anc_chrom_2),
                                  phased = TRUE,
                                  pop_size = 1000,
                                freq_ancestor_1 = 0.1,
                                t = 30)

  ll_900 <- log_likelihood_diploid(cbind(1,
                                        local_data$location,
                                        local_data$anc_chrom_1,
                                        local_data$anc_chrom_2),
                                  phased = TRUE,
                                  pop_size = 1000,
                                  freq_ancestor_1 = 0.1,
                                t = 900)
  testthat::expect_gte(ll_30, ll_900)

  multi_ll <- log_likelihood_diploid(cbind(1,
                                           local_data$location,
                                           local_data$anc_chrom_1,
                                           local_data$anc_chrom_2),
                                     phased = TRUE,
                                     pop_size = 1000,
                                     freq_ancestor_1 = 0.1,
                                   t = c(30, 100, 300))

  testthat::expect_true(length(multi_ll) == 3)
})

test_that("phased, expectation", {
  testthat::skip_on_os("solaris")
  vx <- sim_phased_unphased(pop_size = 10000,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 100,
                            size_in_morgan = 1,
                            time_points = c(20),
                            record_true_junctions = TRUE)

  a <- mean(c(vx$true_results$junctions_chrom_1,
              vx$true_results$junctions_chrom_2))

  expected_j <- number_of_junctions(N = 10000, t = 20)

  testthat::expect_true(abs(a - expected_j) / expected_j < 0.5)
})
