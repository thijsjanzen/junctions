context("unphased")
test_that("unphased, use", {

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


    estim_time <- estimate_time_unphased(cbind(local_data$anc_chrom_1,
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
  ll_100 <- loglikelihood_unphased(cbind(local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   local_data$location,
                                   pop_size = 100,
                                   freq_ancestor_1 = 0.5,
                                   t = 100)

  ll_200 <- loglikelihood_unphased(cbind(local_data$anc_chrom_1,
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

  ll_30 <- loglikelihood_unphased(cbind(local_data$anc_chrom_1,
                                        local_data$anc_chrom_2),
                                  local_data$location,
                                  pop_size = 1000,
                                  freq_ancestor_1 = 0.1,
                                  t = 30)

  ll_100 <- loglikelihood_unphased(cbind(local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   local_data$location,
                                   pop_size = 1000,
                                   freq_ancestor_1 = 0.1,
                                   t = 600)
  testthat::expect_gte(ll_30, ll_100)

  ll_inf <- loglikelihood_unphased(cbind(local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   local_data$location,
                                   pop_size = 1000,
                                   freq_ancestor_1 = 0.1,
                                   t = 0)
  testthat::expect_true(is.infinite(ll_inf))

  multi_ll <- loglikelihood_unphased(cbind(local_data$anc_chrom_1,
                                           local_data$anc_chrom_2),
                                     local_data$location,
                                     pop_size = 1000,
                                     freq_ancestor_1 = 0.1,
                                     t = c(0, 10, 20))
})

test_that("unphased, exceptions", {

  vx <- sim_phased_unphased(pop_size = 100,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 201,
                            size_in_morgan = 1,
                            markers = 1000,
                            seed = 42)

  testthat::expect_warning(
    vx <- sim_phased_unphased(pop_size = 100,
                              freq_ancestor_1 = 0.5,
                              total_runtime = 201,
                              size_in_morgan = 1,
                              markers = 1000,
                              time_points = c(100, 200))
  )
})

test_that("unphased, pop size", {
  vx <- sim_phased_unphased(pop_size = 100,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 1000,
                            size_in_morgan = 1,
                            markers = 1000,
                            seed = 42)

  local_data <- subset(vx, vx$individual == 0 &
                         vx$time == 1000)

  ll_vals <- c()
  n_vals <- seq(0, 500, by = 5)
  for( n in n_vals) {
    local_ll <- loglikelihood_unphased(cbind(local_data$anc_chrom_1,
                                             local_data$anc_chrom_2),
                                       local_data$location,
                                       pop_size = n,
                                       freq_ancestor_1 = 0.1,
                                       t = 1000)
    ll_vals <- c(ll_vals, local_ll)
  }
  a <- n_vals[which.max(ll_vals)]
  plot(exp(ll_vals)~n_vals, type = "l")

  found <- c()
  for( a in unique(vx$individual)) {
    local_data <- subset(vx, vx$individual == a &
                           vx$time == 1000)
    vy <- estimate_time_unphased(cbind(local_data$anc_chrom_1,
                                       local_data$anc_chrom_2),
                                 local_data$location,
                                 freq_ancestor_1 = 0.5,
                                 upper_lim = 2000,
                                 optim_pop_size = TRUE,
                                 verbose = FALSE)
    found <- rbind(found, vy$par)
  }
  estimates <- colMeans(found)
#  testthat::expect_equal(estimates[1], 1000, tolerance = 0.1, scale = 1)
#  testthat::expect_equal(estimates[2], 100, tolerance = 0.1, scale = 1)

}