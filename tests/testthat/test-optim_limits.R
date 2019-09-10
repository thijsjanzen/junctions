context("optim_limits")
test_that("optim_limits, use", {

  vx <- sim_phased_unphased(pop_size = 100,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 101,
                            size_in_morgan = 1,
                            number_of_markers = 1000,
                            time_points = c(100, 200),
                            seed = 42)


  focal_time <- 100
  local_data <- subset(vx, vx$individual == 1 &
                         vx$time == focal_time)

  estim_time <- estimate_time_phased(cbind(local_data$anc_chrom_1,
                                           local_data$anc_chrom_2),
                                     local_data$location,
                                     pop_size = 100,
                                     freq_ancestor_1 = 0.5,
                                     upper_lim = 200)

  estim <- function(lower, upper) {
    estimate_time_phased(cbind(local_data$anc_chrom_1,
                               local_data$anc_chrom_2),
                         local_data$location,
                         pop_size = 100,
                         freq_ancestor_1 = 0.5,
                         lower_lim = lower,
                         upper_lim = upper)
  }

  estim_time2 <- optim_limits(lower = 2, upper = 200, calc_func = estim)

  testthat::expect_equal(estim_time, estim_time2)


  # test scaling up upper limit:
  estim_time3 <- optim_limits(lower = 2, upper = 20, calc_func = estim)
  estim_time3 <- optim_limits(lower = 2, upper = 20, calc_func = estim,
                              iterations = 1)

  # test scaling down lower limit:
  estim_time4 <- optim_limits(lower = 500, upper = 1000, calc_func = estim)


  estim_time <- optim_limits(lower = 500, upper = 1000, calc_func = estim
                             verbose = TRUE)

  estim_time4 <- optim_limits(lower = 500, upper = 1000, calc_func = estim,
                              iterations = 1)



})
