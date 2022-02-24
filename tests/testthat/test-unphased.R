context("unphased")
test_that("unphased, use", {
  testthat::skip_on_os("solaris")
  population_size <- 100
  max_t <- 110
  vx <- sim_phased_unphased(pop_size = population_size,
                            freq_ancestor_1 = 0.5,
                            total_runtime = max_t,
                            size_in_morgan = 1,
                            markers = 1000,
                            time_points = c(50, 100))

  num_indiv <- length(unique(vx$individual))
  testthat::expect_equal(num_indiv, 10)
  testthat::expect_equal(length(unique(vx$time)), 2)

  local_data <- subset(vx, vx$individual == 0 &
                         vx$time == 100)

  ll_100 <- log_likelihood_diploid(cbind(1,
                                         local_data$location,
                                         local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   pop_size = 100,
                                   freq_ancestor_1 = 0.5,
                                   t = 100,
                                   phased = FALSE)

  ll_2000 <- log_likelihood_diploid(cbind(1,
                                         local_data$location,
                                         local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   pop_size = 100,
                                   freq_ancestor_1 = 0.5,
                                   phased = FALSE,
                                   t = 2000)

  testthat::expect_gte(ll_100, ll_2000)

  vx <- sim_phased_unphased(pop_size = 10000,
                            freq_ancestor_1 = 0.1,
                            total_runtime = 20,
                            size_in_morgan = 1,
                            markers = 10000,
                            time_points = c(20))

  local_data <- subset(vx, vx$individual == 0 &
                         vx$time == 20)

  ll_30 <- log_likelihood_diploid(cbind(1, local_data$location,
                                        local_data$anc_chrom_1,
                                        local_data$anc_chrom_2),
                                  pop_size = 1000,
                                  freq_ancestor_1 = 0.1,
                                  phased = FALSE,
                                  t = 30)

  ll_100 <- log_likelihood_diploid(cbind(1, local_data$location,
                                         local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   pop_size = 1000,
                                   freq_ancestor_1 = 0.1,
                                   phased = FALSE,
                                   t = 600)
  testthat::expect_gte(ll_30, ll_100)

  multi_ll <- log_likelihood_diploid(cbind(1,
                                           local_data$location,
                                           local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   pop_size = 1000,
                                   freq_ancestor_1 = 0.1,
                                   phased = FALSE,
                                   t = c(3, 30000))

  testthat::expect_true(length(multi_ll) == 2)
  testthat::expect_gt(multi_ll[1], multi_ll[2])
})

test_that("unphased, time points", {
  testthat::skip_on_os("solaris")
  population_size <- 100
  max_t <- 10
  vx <- sim_phased_unphased(pop_size = population_size,
                            freq_ancestor_1 = 0.5,
                            total_runtime = max_t,
                            size_in_morgan = 1,
                            markers = 1000,
                            time_points = -1)

  sim_t <- unique(vx$time)
  testthat::expect_equal(length(sim_t), max_t + 1)  # [0, 1, ..., max_t]

  testthat::expect_warning(
    vx <- sim_phased_unphased(pop_size = population_size,
                              freq_ancestor_1 = 0.5,
                              total_runtime = max_t,
                              size_in_morgan = 1,
                              markers = 1000,
                              time_points = max_t + 5)

  )
  testthat::expect_equal(length(unique(vx$time)), 1)
})

test_that("unphased, junctions", {
  testthat::skip_on_os("solaris")
  N <- 10000 # nolint
  R <- 10000 # nolint
  t <- 10
  H_0 <- 0.5 # nolint
  C <- 1     # nolint

  testthat::expect_output(
    vx <- sim_phased_unphased(pop_size = N,
                              freq_ancestor_1 = H_0,
                              total_runtime = t,
                              size_in_morgan = C,
                              markers = R,
                              time_points = t,
                              num_indiv_sampled = 100,
                              record_true_junctions = TRUE,
                              verbose = TRUE)
  )

  num_j_true <- mean(c(vx$true_results$junctions_chrom_1,
                       vx$true_results$junctions_chrom_2))
  vx <- vx$results

  num_j <- c()
  for (i in unique(vx$individual)) {
    dd <- subset(vx, vx$individual == i)
    num_j_1 <- sum(abs(diff(dd$anc_chrom_1)))
    num_j_2 <- sum(abs(diff(dd$anc_chrom_1)))
    num_j <- c(num_j, c(num_j_1, num_j_2))
  }

  obs_j <- mean(num_j)
  exp_j <- junctions::number_of_junctions(N = N,
                                          R = R,
                                          H_0 = H_0,
                                          C = C,
                                          t = t)

  testthat::expect_equal(obs_j, exp_j, tolerance = 0.2)

  N <- 10000 # nolint
  R <- 10000 # nolint
  t <- 20
  H_0 <- 0.5 # nolint
  C <- 1    # nolint

  vx <- sim_phased_unphased(pop_size = N,
                            freq_ancestor_1 = H_0,
                            total_runtime = t,
                            size_in_morgan = C,
                            markers = R,
                            time_points = t,
                            num_indiv_sampled = 30)

  num_j <- c()
  for (i in unique(vx$individual)) {
    dd <- subset(vx, vx$individual == i)
    num_j_1 <- sum(abs(diff(dd$anc_chrom_1)))
    num_j_2 <- sum(abs(diff(dd$anc_chrom_1)))
    num_j <- c(num_j, c(num_j_1, num_j_2))
  }

  obs_j <- mean(num_j)
  exp_j <- junctions::number_of_junctions(N = N,
                                          R = R,
                                          H_0 = H_0,
                                          C = C,
                                          t = t)

  testthat::expect_equal(obs_j, exp_j, tolerance = 0.2)
})
