context("phased")
test_that("phased, use", {

  vx <- sim_phased_unphased(pop_size = 100,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 201,
                            size_in_morgan = 1,
                            markers = 1000,
                            time_points = c(100, 200),
                            seed = 42)

  num_indiv <- length(unique(vx$individual))
  testthat::expect_equal(num_indiv, 10)
  testthat::expect_equal(length(unique(vx$time)), 2)

  found <- c()
  focal_time <- 100
  for (indiv in 0:9) {
    local_data <- subset(vx, vx$individual == indiv &
                           vx$time == focal_time)


    estim_time <- estimate_time_phased(cbind(local_data$anc_chrom_1,
                                             local_data$anc_chrom_2),
                                       local_data$location,
                                       pop_size = 100,
                                       freq_ancestor_1 = 0.5,
                                       upper_lim = 200)
    found[indiv + 1] <- estim_time$minimum
  }

  testthat::expect_equal(mean(found), focal_time, tolerance = 10)



  local_data <- subset(vx, vx$individual == 0 &
                         vx$time == 100)
  ll_100 <- loglikelihood_phased(cbind(local_data$anc_chrom_1,
                                       local_data$anc_chrom_2),
                                 local_data$location,
                                 pop_size = 100,
                                 freq_ancestor_1 = 0.5,
                                 t = 100)

  ll_200 <- loglikelihood_phased(cbind(local_data$anc_chrom_1,
                                       local_data$anc_chrom_2),
                                 local_data$location,
                                 pop_size = 100,
                                 freq_ancestor_1 = 0.5,
                                 t = 200)

  testthat::expect_gte(ll_100, ll_200)

  vx <- sim_phased_unphased(pop_size = 1000,
                            freq_ancestor_1 = 0.1,
                            total_runtime = 30,
                            size_in_morgan = 1,
                            markers = 1000,
                            time_points = c(30),
                            seed = 421)

  local_data <- subset(vx, vx$individual == 0 &
                         vx$time == 30)

  ll_30 <- loglikelihood_phased(cbind(local_data$anc_chrom_1,
                                      local_data$anc_chrom_2),
                                local_data$location,
                                pop_size = 1000,
                                freq_ancestor_1 = 0.1,
                                t = 30)

  ll_90 <- loglikelihood_phased(cbind(local_data$anc_chrom_1,
                                      local_data$anc_chrom_2),
                                local_data$location,
                                pop_size = 1000,
                                freq_ancestor_1 = 0.1,
                                t = 90)
  testthat::expect_gte(ll_30, ll_90)

  ll_inf <- loglikelihood_phased(cbind(local_data$anc_chrom_1,
                                       local_data$anc_chrom_2),
                                 local_data$location,
                                 pop_size = 1000,
                                 freq_ancestor_1 = 0.1,
                                 t = 0)
  testthat::expect_true(is.infinite(ll_inf))


  multi_ll <- loglikelihood_phased(cbind(local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   local_data$location,
                                   pop_size = 1000,
                                   freq_ancestor_1 = 0.1,
                                   t = c(0, 10, 20))
})

test_that("phased, pop size", {
  vx <- sim_phased_unphased(pop_size = 100,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 501,
                            size_in_morgan = 1,
                            markers = 1000,
                            time_points = c(500),
                            seed = 42)

  focal_data <- subset(vx, vx$individual == 0 & vx$time == 500)
  phased_data <- cbind(focal_data$anc_chrom_1, focal_data$anc_chrom_2)
  time2 <- estimate_time_unphased(local_anc_matrix = phased_data,
                                  locations = focal_data$location,
                                  pop_size = 500,
                                  freq_ancestor_1 = 0.5,
                                  optim_pop_size = TRUE,
                                  verbose = TRUE)

  testthat::expect_true(abs(time2$par[1] - 500) / 500 < 0.5)
  testthat::expect_true(abs(time2$par[2] - 100) / 100 < 0.5)

  time3 <- estimate_time_phased(local_anc_matrix = phased_data,
                                locations = focal_data$location,
                                pop_size = 500,
                                freq_ancestor_1 = 0.5,
                                optim_pop_size = TRUE,
                                verbose = TRUE)

  testthat::expect_true(abs(time3$par[1] - 500) / 500 < 0.5)
  testthat::expect_true(abs(time3$par[2] - 100) / 100 < 0.5)

  testthat::expect_false(time2$par[1] == time3$par[1])
})


test_that("phased, abuse", {
  vx <- sim_phased_unphased(pop_size = 100,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 101,
                            size_in_morgan = 1,
                            markers = 1000,
                            time_points = c(100),
                            seed = 42)

  focal_data <- subset(vx, vx$individual == 0 & vx$time == 100)

  alt_data <- rep(NA, length(focal_data$location))
  alt_data[focal_data$anc_chrom_1 == 0] <- 1
  alt_data[focal_data$anc_chrom_2 == 1] <- 2
  alt_data[focal_data$anc_chrom_1 != focal_data$anc_chrom_2] <- 3

  time1 <- estimate_time_unphased(local_anc_matrix = alt_data,
                                  locations = focal_data$location,
                                  pop_size = 100,
                                  freq_ancestor_1 = 0.5,
                                  verbose = TRUE)

  phased_data <- cbind(focal_data$anc_chrom_1, focal_data$anc_chrom_2)
  time2 <- estimate_time_unphased(local_anc_matrix = phased_data,
                                  locations = focal_data$location,
                                  pop_size = 100,
                                  freq_ancestor_1 = 0.5)

  testthat::expect_equal(time1, time2)

  alt_data[alt_data == 3] <- 5
  testthat::expect_error(
    estimate_time_phased(local_anc_matrix = alt_data,
                         locations = focal_data$location,
                         pop_size = 100,
                         freq_ancestor_1 = 0.5)
  )

  testthat::expect_error(
    loglikelihood_unphased(local_anc_matrix = alt_data,
                           locations = focal_data$location,
                           pop_size = 100,
                           freq_ancestor_1 = 0.5,
                           t = 10)
  )

  testthat::expect_error(
    loglikelihood_phased(local_anc_matrix = alt_data,
                         locations = focal_data$location,
                         pop_size = 100,
                         freq_ancestor_1 = 0.5,
                         t = 10)
  )



  alt_data <- rep(NA, length(focal_data$location))
  alt_data[focal_data$anc_chrom_1 == 0] <- 1
  alt_data[focal_data$anc_chrom_2 == 1] <- 2
  alt_data[focal_data$anc_chrom_1 == 0 & focal_data$anc_chrom_2 == 1] <- 3
  alt_data[focal_data$anc_chrom_1 == 1 & focal_data$anc_chrom_2 == 0] <- 4

  testthat::expect_error(estimate_time_unphased(local_anc_matrix = alt_data,
                                                locations = focal_data$location,
                                                pop_size = 100,
                                                freq_ancestor_1 = 0.5))

  testthat::expect_silent(estimate_time_phased(local_anc_matrix = alt_data,
                                               locations = focal_data$location,
                                               pop_size = 100,
                                               freq_ancestor_1 = 0.5))
})

test_that("phased, expectation", {

  vx <- sim_phased_unphased(pop_size = 10000,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 100,
                            size_in_morgan = 1,
                            time_points = c(20),
                            seed = 42,
                            record_true_junctions = TRUE)
  a <- mean(c(vx$true_results$junctions_chrom_1,
              vx$true_results$junctions_chrom_2))

  expected_j <- number_of_junctions(N = 10000, t = 20)

  testthat::expect_true( abs(a - expected_j) / expected_j < 0.1)
})
