context("estimate_time_haploid")
test_that("estimate_time_haploid, use", {
  testthat::skip_on_os("solaris")
  sim_time <- 100
  vx <- sim_phased_unphased(total_runtime = sim_time,
                            time_points = sim_time,
                            pop_size = 1000,
                            num_indiv_sampled = 20)

  estim_time <- estimate_time_haploid(cbind(vx$individual,
                                            vx$location,
                                            vx$anc_chrom_1),
                                      N = 1000,
                                      freq_ancestor_1 = 0.5)

  testthat::expect_length(estim_time$time, 1)
  testthat::expect_length(estim_time$loglikelihood, 1)


  ll1 <- log_likelihood_haploid(ancestry_matrix = cbind(vx$individual,
                                                        vx$location,
                                                        vx$anc_chrom_1),
                                N = 1000,
                                freq_ancestor_1 = 0.5,
                                t = 100)

  testthat::expect_length(ll1, 1)

  focal_t <- 85:100
  ll <- log_likelihood_haploid(ancestry_matrix = cbind(vx$individual,
                                                       vx$location,
                                                       vx$anc_chrom_1),
                               N = 1000,
                               freq_ancestor_1 = 0.5,
                               t = focal_t)

  testthat::expect_length(ll, length(focal_t))

  # check boundaries:
  testthat::expect_warning(
    estimate_time_haploid(cbind(vx$individual,
                                vx$location,
                                vx$anc_chrom_1),
                          N = 1000,
                          freq_ancestor_1 = 0.5,
                          upper_lim = 10)
  )

  testthat::expect_output(
    estimate_time_haploid(cbind(vx$individual,
                                vx$location,
                                vx$anc_chrom_1),
                          N = 1000,
                          freq_ancestor_1 = 0.5,
                          verbose = TRUE)
  )

  # check it works with tibble:
  testthat::expect_silent(
    estimate_time_haploid(tibble::tibble(vx$individual,
                                vx$location,
                                vx$anc_chrom_1),
                          N = 1000,
                          freq_ancestor_1 = 0.5)
  )

})

test_that("estimate_time_diploid, use", {
  testthat::skip_on_os("solaris")
  sim_time <- 100
  vx <- sim_phased_unphased(total_runtime = sim_time,
                            time_points = sim_time,
                            pop_size = 1000,
                            num_indiv_sampled = 5)

  indiv <- vx$individual
  indiv[indiv < 2] <- 1
  indiv[indiv >= 2] <- 2

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
                              analysis_type = "separate",
                              pop_size = 1000,
                              freq_ancestor_1 = 0.5)

  # chromosomes and "all" is equal in this case,
  # because we did not add multiple individuals
  testthat::expect_equal(t2$loglikelihood, t3$loglikelihood)
  testthat::expect_equal(t2$time, t3$time)

  # 30 chromosomes:
  testthat::expect_equal(length(t2$time), 5)
  testthat::expect_equal(length(t3$time), 5)
  # 2 individuals:
  testthat::expect_equal(length(t1$time), 2)


  # test verbose output:
  testthat::expect_output(
    estimate_time_diploid(ancestry_information = cbind(indiv, vx$individual,
                                                             vx$location,
                                                             vx$anc_chrom_1,
                                                             vx$anc_chrom_2),
                                analysis_type = "individuals",
                                pop_size = 1000,
                                freq_ancestor_1 = 0.5,
                                verbose = TRUE)
  )

  # test error:
  testthat::expect_error(
    estimate_time_diploid(ancestry_information = cbind(indiv, vx$individual,
                                                       vx$location,
                                                       vx$anc_chrom_1,
                                                       vx$anc_chrom_2),
                          analysis_type = "individual",
                          pop_size = 1000,
                          freq_ancestor_1 = 0.5)
  )

  # test working with tibble
  input_tibble <- tibble::tibble(indiv, vx$individual,
                         vx$location,
                         vx$anc_chrom_1,
                         vx$anc_chrom_2)
  testthat::expect_silent(
    estimate_time_diploid(ancestry_information = input_tibble,
                          analysis_type = "individuals",
                          pop_size = 1000,
                          freq_ancestor_1 = 0.5)
  )
})
