context("regressions")

test_that("simulation input is validated instead of crashing", {
  testthat::expect_error(sim_phased_unphased(total_runtime = 5,
                                             time_points = -5),
                         "non-negative")
  testthat::expect_error(sim_phased_unphased(total_runtime = 5,
                                             time_points = 2.5),
                         "whole numbers")
  testthat::expect_error(sim_phased_unphased(pop_size = 1), "two individuals")
  testthat::expect_error(sim_phased_unphased(pop_size = 10,
                                             num_indiv_sampled = 20),
                         "exceed population size")
  testthat::expect_error(sim_backcrossing(number_of_markers = 0))
  testthat::expect_error(sim_backcrossing(freq_ancestor_1 = 0.5),
                         "deprecated")
  # used to hang: more crossovers than recombination sites
  vx <- sim_fin_chrom(R = 1, morgan = 1, total_runtime = 5)
  testthat::expect_equal(length(vx$avgJunctions), 5)
})

test_that("num_threads = -1 does not hang", {
  vx <- sim_phased_unphased(total_runtime = 3, num_threads = -1, seed = 1)
  testthat::expect_gt(nrow(vx), 0)
})

test_that("sim_phased_unphased is reproducible", {
  a <- sim_phased_unphased(total_runtime = 10, time_points = 10, seed = 5)
  b <- sim_phased_unphased(total_runtime = 10, time_points = 10, seed = 5,
                           num_threads = 2)
  testthat::expect_identical(a, b)

  set.seed(3)
  c1 <- sim_phased_unphased(total_runtime = 10, time_points = 10)
  set.seed(3)
  c2 <- sim_phased_unphased(total_runtime = 10, time_points = 10)
  testthat::expect_identical(c1, c2)

  # a given seed leaves the global RNG state untouched
  set.seed(9)
  r1 <- stats::runif(1)
  set.seed(9)
  sim_phased_unphased(total_runtime = 3, time_points = 3, seed = 1)
  r2 <- stats::runif(1)
  testthat::expect_equal(r1, r2)

  # the fin/inf simulations no longer reset the global RNG
  set.seed(9)
  r1 <- stats::runif(1)
  set.seed(9)
  sim_inf_chrom(total_runtime = 3)
  r2 <- stats::runif(1)
  testthat::expect_equal(r1, r2)
})

test_that("user supplied markers are used as given", {
  vx <- sim_phased_unphased(pop_size = 200, total_runtime = 20,
                            time_points = 20,
                            markers = c(0.6, 0.2, 1.0),
                            num_indiv_sampled = 200, seed = 2)
  testthat::expect_equal(sort(unique(vx$location)), c(0.2, 0.6, 1.0))
  # a marker at the very end of the chromosome used to always read as 0
  freq <- tapply(vx$anc_chrom_1, vx$location, mean)
  testthat::expect_gt(freq[["1"]], 0.3)
})

test_that("phasing error only flips the expected number of markers", {
  vx <- sim_phased_unphased(pop_size = 100, total_runtime = 5,
                            time_points = 5, markers = 1,
                            num_indiv_sampled = 100,
                            error_rate = 1e-6, seed = 1)
  testthat::expect_equal(vx$true_data$anc_chrom_1[vx$true_data$location %in%
                                                    vx$phased_data$location],
                         vx$phased_data$anc_chrom_1)
})

test_that("time estimation handles edge cases", {
  vx <- sim_phased_unphased(pop_size = 1000, total_runtime = 50,
                            markers = 200, time_points = 50,
                            num_indiv_sampled = 1, seed = 1)
  d <- cbind(1, 1, vx$location, vx$anc_chrom_1, vx$anc_chrom_2)

  # bounds that do not contain the old fixed starting value
  t1 <- estimate_time_diploid(d, lower_lim = 20, upper_lim = 200)
  testthat::expect_true(is.finite(t1$loglikelihood))
  testthat::expect_gte(t1$time, 20)
  testthat::expect_error(estimate_time_diploid(d, lower_lim = 50,
                                               upper_lim = 10))

  # single marker chromosomes are dropped with a warning, not a crash
  d2 <- rbind(d, c(1, 2, 0.5, 0, 1))
  testthat::expect_warning(t2 <- estimate_time_diploid(d2,
                                                       analysis_type = "all"),
                           "single marker")
  testthat::expect_equal(t2$loglikelihood,
                         estimate_time_diploid(d, analysis_type = "all")$
                           loglikelihood)
  testthat::expect_error(estimate_time_diploid(d[1, , drop = FALSE]),
                         "at least two markers")

  # ancestry has to be 0/1
  d3 <- d
  d3[, 4:5] <- d3[, 4:5] + 1
  testthat::expect_error(estimate_time_diploid(d3), "coded as 0")
  d3 <- d
  d3[1, 4] <- NA
  testthat::expect_error(estimate_time_diploid(d3), "NA")

  # character identifiers are kept in the output
  df <- data.frame(id = "fishA", chr = "chr1", pos = d[, 3],
                   a = d[, 4], b = d[, 5])
  t4 <- estimate_time_diploid(df, analysis_type = "separate")
  testthat::expect_equal(t4$individual, "fishA")
  testthat::expect_equal(t4$chromosome, "chr1")
  testthat::expect_equal(t4$time, estimate_time_diploid(d)$time)

  # "all": chromosomes of different individuals are never merged, even if
  # chromosome labels and positions would allow it
  d5 <- rbind(cbind(1, 1, d[, 3] / 2, d[, 4:5]),
              cbind(2, 1, 0.5 + d[, 3] / 2, d[, 4:5]))
  t5 <- estimate_time_diploid(d5, analysis_type = "all")
  t6 <- estimate_time_diploid(d5, analysis_type = "individuals")
  testthat::expect_equal(t5$loglikelihood, sum(t6$loglikelihood),
                         tolerance = 1e-6)

  # haploid: single marker chromosome
  hm <- cbind(c(1, 1, 1, 2), c(0.1, 0.2, 0.3, 0.5), c(0, 1, 1, 0))
  testthat::expect_silent(log_likelihood_haploid(hm, N = 100, t = 10))
})

test_that("C++ errors reach R with their message", {
  m <- cbind(1, c(0.1, 0.2, 0.3), c(0, 0, 1), c(0, 1, 1))
  testthat::expect_error(log_likelihood_diploid(m, pop_size = 100, t = 0.5),
                         "t < 1")
  testthat::expect_error(log_likelihood_diploid(m, pop_size = 1, t = 10),
                         "pop_size < 2")
  # thrown inside the nlopt objective
  testthat::expect_error(estimate_time_diploid(cbind(1, m),
                                               freq_ancestor_1 = 1),
                         "p >= 1")
})

test_that("analytical functions are vectorised", {
  n <- c(100, 1000, Inf)
  r <- c(1000, Inf, Inf)
  testthat::expect_equal(number_of_junctions(N = n, R = r, t = 10),
                         vapply(1:3, function(i)
                           number_of_junctions(N = n[i], R = r[i], t = 10),
                           numeric(1)))
  testthat::expect_equal(calc_k(N = c(Inf, 100), R = c(100, Inf)),
                         c(50, 100))
  testthat::expect_warning(et <- estimate_time(J = c(5, 500), N = 100))
  testthat::expect_true(is.na(et[2]))
  testthat::expect_equal(et[1], estimate_time(J = 5, N = 100))
  testthat::expect_warning(te <- time_error(t = c(10, 1e5), N = 100,
                                            R = 1000))
  testthat::expect_true(is.na(te[2]))
  testthat::expect_equal(length(calculate_mat(N = c(100, 1000), R = 1000)), 2)
  testthat::expect_equal(number_of_junctions_markers(N = 100, t = 50,
                           marker_distribution = c(0.9, 0.1, 0.5)),
                         number_of_junctions_markers(N = 100, t = 50,
                           marker_distribution = c(0.1, 0.5, 0.9)))
})

test_that("diploid log likelihood matches reference values", {
  # reference values computed before the task_arena change; guards against
  # e.g. the first (unconditional) marker pair being overwritten
  m <- cbind(c(1, 1, 1, 1, 1, 2, 2, 2),
             c(0.05, 0.1, 0.3, 0.32, 0.8, 0.1, 0.4, 0.9),
             c(0, 0, 1, 1, 0, 1, 1, 0),
             c(0, 1, 1, 1, 1, 0, 1, 0))
  ll <- function(t, phased, threads = 1) {
    log_likelihood_diploid(m, pop_size = 1000, t = t, phased = phased,
                           num_threads = threads)
  }
  testthat::expect_equal(ll(5, FALSE), -8.7759337170, tolerance = 1e-8)
  testthat::expect_equal(ll(50, FALSE), -8.3554718334, tolerance = 1e-8)
  testthat::expect_equal(ll(5, TRUE), -10.8553752587, tolerance = 1e-8)
  testthat::expect_equal(ll(50, TRUE), -10.4349133750, tolerance = 1e-8)
  # same answer regardless of the number of threads
  testthat::expect_equal(ll(50, FALSE, threads = 2), ll(50, FALSE))
  testthat::expect_equal(ll(50, FALSE, threads = -1), ll(50, FALSE))
})

test_that("chromosomes stacked under one label are split on location", {
  # as in the vignette: several individuals passed with the same individual
  # and chromosome label; a decrease in location starts a new chromosome
  vx <- sim_phased_unphased(pop_size = 1000, total_runtime = 50,
                            markers = 200, time_points = 50,
                            num_indiv_sampled = 3, seed = 3)
  stacked <- cbind(1, 1, vx$location, vx$anc_chrom_1, vx$anc_chrom_2)
  labelled <- cbind(vx$individual, 1, vx$location,
                    vx$anc_chrom_1, vx$anc_chrom_2)

  t_stacked <- estimate_time_diploid(stacked, analysis_type = "all")
  t_labelled <- estimate_time_diploid(labelled, analysis_type = "all")
  testthat::expect_true(is.finite(t_stacked$loglikelihood))
  testthat::expect_equal(t_stacked, t_labelled)

  testthat::expect_equal(
    log_likelihood_diploid(stacked[, 2:5], pop_size = 1000, t = 50),
    sum(vapply(unique(vx$individual), function(i) {
      log_likelihood_diploid(labelled[labelled[, 1] == i, 2:5],
                             pop_size = 1000, t = 50)
    }, numeric(1))))
})
