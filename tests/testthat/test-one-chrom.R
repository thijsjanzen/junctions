context("one_chromosome")
test_that("one chrom, use", {

  population_size <- 100000
  run_time <- 10

  vx <- sim_phased_unphased(pop_size = population_size,
                            total_runtime = run_time,
                            markers = 10000,
                            seed = 4,
                            time_points = run_time)

  focal_data <- subset(vx, vx$time == run_time & vx$individual == 0)
  time1 <- estimate_time_one_chrom(J = sum(abs(diff(focal_data$anc_chrom_1))),
                                   N = population_size,
                                   H_0 = 0.5,
                                   marker_distribution = focal_data$location)
  time2 <- estimate_time_one_chrom(J = sum(abs(diff(focal_data$anc_chrom_2))),
                                   N = population_size,
                                   H_0 = 0.5,
                                   marker_distribution = focal_data$location)
  cat(time1, time2, "\n")
  testthat::expect_true( abs((time1 + time2) / 2 - run_time) < 5)

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