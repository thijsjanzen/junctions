context("one_chromosome")
test_that("one chrom, use", {
  testthat::skip_on_os("solaris")
  population_size <- 10000
  run_time <- 10

  vx <- sim_phased_unphased(pop_size = population_size,
                            total_runtime = run_time,
                            markers = 1000,
                            time_points = run_time,
                            num_indiv_sampled = 5)

  found <- c()

  for (i in unique(vx$individual)) {
    focal_data <- subset(vx, vx$time == run_time & vx$individual == i)
    time1 <- estimate_time_one_chrom(J = sum(abs(diff(focal_data$anc_chrom_1))),
                                     N = population_size,
                                     H_0 = 0.5,
                                     marker_distribution = focal_data$location)
    time2 <- estimate_time_one_chrom(J = sum(abs(diff(focal_data$anc_chrom_2))),
                                     N = population_size,
                                     H_0 = 0.5,
                                     marker_distribution = focal_data$location)

    found <- c(found, c(time1, time2))
  }

  testthat::expect_equal(mean(found), run_time, tolerance = 3)

  # induce marker error
  testthat::expect_error(
     estimate_time_one_chrom(J = sum(abs(diff(focal_data$anc_chrom_1))),
                                   N = 1000,
                                   H_0 = 0.5)
  )

  # induce warning:
  population_size <- 100
  run_time <- 100

  vx <- sim_phased_unphased(pop_size = population_size,
                            total_runtime = run_time,
                            markers = 1000,
                            time_points = run_time,
                            num_indiv_sampled = 2)
  focal_data <- subset(vx, vx$time == run_time & vx$individual == 0)
  testthat::expect_warning(
    estimate_time_one_chrom(J = sum(abs(diff(focal_data$anc_chrom_2))),
                            N = population_size,
                            H_0 = 0.5,
                            marker_distribution = focal_data$location,
                            upper_lim = 10)
  )
})
