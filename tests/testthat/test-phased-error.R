context("phased_with_error")
test_that("phased, use", {
  vx <- sim_phased_unphased(pop_size = 100,
                              freq_ancestor_1 = 0.5,
                              total_runtime = 21,
                              size_in_morgan = 1,
                              markers = 1000,
                              time_points = c(10, 20),
                              coverage = 0.99,
                              error_rate = 0)

  num_indiv <- length(unique(vx$true_data$individual))
  testthat::expect_equal(num_indiv, 10)
  testthat::expect_equal(length(unique(vx$true_data$time)), 2)

  # now we introduce less coverage:
  vx <- sim_phased_unphased(pop_size = 100,
                              freq_ancestor_1 = 0.5,
                              total_runtime = 21,
                              size_in_morgan = 1,
                              markers = 1000,
                              time_points = c(10, 20),
                              coverage = 0.5,
                              error_rate = 0)

  true_data <- vx$true_data
  phased_data <- vx$phased_data
  true_markers <- length(unique(true_data$location))
  phased_markers <- length(unique(phased_data$location))
  testthat::expect_true(phased_markers / true_markers == 0.5)

  # now we introduce error
  errorrr <- 0.5
  vx <- sim_phased_unphased(pop_size = 10000,
                              freq_ancestor_1 = 0.5,
                              total_runtime = 20,
                              size_in_morgan = 1,
                              markers = 1000,
                              time_points = c(20),
                              coverage = 1,
                              error_rate = errorrr)

  true_data <- vx$true_data
  phased_data <- vx$phased_data
  expected_heterozygosity <- 2 * 0.5 * 0.5 * (1 - 1 / (2 * 100)) ^ 200

  for (i in unique(true_data$individual)) {
    a <- subset(true_data, true_data$individual == i)
    b <- subset(phased_data, phased_data$individual == i)

    a1 <- sum(a$anc_chrom_1 != b$anc_chrom_1)
    a2 <- sum(a$anc_chrom_2 != b$anc_chrom_2)

    testthat::expect_equal(a1,
                           errorrr *
                             expected_heterozygosity * length(a$location),
                           tolerance = 10)
    testthat::expect_equal(a2,
                           errorrr *
                             expected_heterozygosity * length(a$location),
                           tolerance = 10)
  }
})
