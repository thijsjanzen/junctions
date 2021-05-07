context("explicit")
test_that("explicit, use", {

  marker_vec <- sort(runif(10000, 0, 1))

  max_gen <- 20

  s1 <- sim_phased_unphased(pop_size = 1000,
                            freq_ancestor_1 = 0.5,
                            total_runtime = max_gen,
                            markers = marker_vec,
                            time_points = max_gen,
                            num_indiv_sampled = 100,
                            use_explicit = FALSE)

  s2 <- sim_phased_unphased(pop_size = 1000,
                            freq_ancestor_1 = 0.5,
                            total_runtime = max_gen,
                            markers = marker_vec,
                            time_points = max_gen,
                            num_indiv_sampled = 100,
                            use_explicit = TRUE)

  testthat::expect_equal(s1$location, s2$location)

  num_j_s1 <- c()
  num_j_s2 <- c()
  for (i in unique(s1$individual)) {

    a <- subset(s1, s1$individual == i)
    j1 <- sum(abs(diff(a$anc_chrom_1)))
    j2 <- sum(abs(diff(a$anc_chrom_2)))
    num_j_s1 <- c(num_j_s1, j1, j2)

    b <- subset(s2, s2$individual == i)
    j1 <- sum(abs(diff(b$anc_chrom_1)))
    j2 <- sum(abs(diff(b$anc_chrom_2)))
    num_j_s2 <- c(num_j_s2, j1, j2)
  }
  pp <- t.test(num_j_s1, num_j_s2)

  testthat::expect_gt(pp$p.value, 0.05)
})
