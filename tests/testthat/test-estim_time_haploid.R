context("estimate_time_haploid")
test_that("estimate_time_haploid, use", {
  sim_time <- 100
  vx <- sim_phased_unphased(total_runtime = sim_time,
                            time_points = sim_time,
                            pop_size = 1000,
                            seed = 1,
                            num_indiv_sampled = 200)

  estim_time <- estimate_time_haploid(cbind(vx$individual,
                                            vx$location,
                                            vx$anc_chrom_1),
                                      N = 1000,
                                      freq_ancestor_1 = 0.5)

  testthat::expect_equal(estim_time$time, sim_time, tolerance = 0.1)
})

test_that("estimate_time_diploid, use", {
  sim_time <- 100
  vx <- sim_phased_unphased(total_runtime = sim_time,
                            time_points = sim_time,
                            pop_size = 1000,
                            seed = 1,
                            num_indiv_sampled = 30)

  indiv <- vx$individual
  indiv[indiv < 15] <- 1
  indiv[indiv >= 15] <- 2

  t1 <- estimate_time_diploid(ancestry_information = cbind(indiv, vx$individual,
                                                           vx$location,
                                                           vx$anc_chrom_1,
                                                           vx$anc_chrom_2),
                              analysis_type = "individuals",
                              pop_size = 1000,
                              freq_ancestor_1 = 0.5)

  t2 <- estimate_time_diploid(ancestry_information = cbind(indiv, vx$individual,
                                                           vx$location,
                                                           vx$anc_chrom_1,
                                                           vx$anc_chrom_2),
                              analysis_type = "chromosomes",
                              pop_size = 1000,
                              freq_ancestor_1 = 0.5)

  t3 <- estimate_time_diploid(ancestry_information = cbind(indiv,
                                                           vx$individual,
                                                           vx$location,
                                                           vx$anc_chrom_1,
                                                           vx$anc_chrom_2),
                              analysis_type = "all",
                              pop_size = 1000,
                              freq_ancestor_1 = 0.5)

  # chromosomes and "all" is equal in this case,
  # because we did not add multiple individuals
  testthat::expect_equal(t2$loglikelihood, t3$loglikelihood)
  testthat::expect_equal(t2$time, t3$time)

  # 30 chromosomes:
  testthat::expect_equal(length(t2$time), 30)
  testthat::expect_equal(length(t3$time), 30)
  # 2 individuals:
  testthat::expect_equal(length(t1$time), 2)
})