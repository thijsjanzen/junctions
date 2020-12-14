context("one_chromosome")
test_that("one chrom, use", {

  vx <- sim_phased_unphased(pop_size = 1000,
                            total_runtime = 100,
                            markers = 1000,
                            seed = 42,
                            time_points = 100)
  focal_data <- subset(vx, vx$time == 100 & vx$individual == 0)
  time1 <- estimate_time_one_chrom(J = sum(abs(diff(focal_data$anc_chrom_1))),
                                   N = 1000,
                                   H_0 = 0.5,
                                   marker_distribution = focal_data$location)
  time2 <- estimate_time_one_chrom(J = sum(abs(diff(focal_data$anc_chrom_2))),
                                   N = 1000,
                                   H_0 = 0.5,
                                   marker_distribution = focal_data$location)

  testthat::expect_true( (time1 + time2) / 2 - 100 < 10)

  # induce marker error
  testthat::expect_error(
     estimate_time_one_chrom(J = sum(abs(diff(focal_data$anc_chrom_1))),
                                   N = 1000,
                                   H_0 = 0.5)
  )

  testthat::expect_warning(
    estimate_time_one_chrom(J = sum(abs(diff(focal_data$anc_chrom_1))),
                            N = 1000,
                            H_0 = 0.5,
                            marker_distribution = focal_data$location,
                            upper_lim = 10)
  )
})