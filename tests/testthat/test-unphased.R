context("unphased")
test_that("unphased, use", {

  population_size <- 100
  max_t <- 100
  vx <- sim_phased_unphased(pop_size = population_size,
                            freq_ancestor_1 = 0.5,
                            total_runtime = max_t,
                            size_in_morgan = 1,
                            markers = 1000,
                            time_points = c(50, 100),
                            seed = 43)

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

  ll_200 <- log_likelihood_diploid(cbind(1,
                                         local_data$location,
                                         local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   pop_size = 100,
                                   freq_ancestor_1 = 0.5,
                                   phased = FALSE,
                                   t = 200)

  testthat::expect_gte(ll_100, ll_200)

  vx <- sim_phased_unphased(pop_size = 1000,
                            freq_ancestor_1 = 0.1,
                            total_runtime = 30,
                            size_in_morgan = 1,
                            markers = 10000,
                            time_points = c(30),
                            seed = 421)

  local_data <- subset(vx, vx$individual == 0 &
                         vx$time == 30)

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
                                   t = c(3, 30000, 300000))

  testthat::expect_true(length(multi_ll) == 3)
  testthat::expect_gt(multi_ll[1], multi_ll[2])
  testthat::expect_gt(multi_ll[2], multi_ll[3])
})

test_that("unphased, exceptions", {

  vx <- sim_phased_unphased(pop_size = 100,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 201,
                            size_in_morgan = 1,
                            markers = 1000,
                            seed = 42)

  testthat::expect_warning(
    vx <- sim_phased_unphased(pop_size = 100,
                              freq_ancestor_1 = 0.5,
                              total_runtime = 201,
                              size_in_morgan = 1,
                              markers = 1000,
                              time_points = c(100, 200))
  )
})

test_that("unphased, junctions", {

  N <- 1000
  R <- 1000
  t <- 100
  H_0 <- 0.5
  C <- 1

  vx <- sim_phased_unphased(pop_size = N,
                            freq_ancestor_1 = H_0,
                            total_runtime = t,
                            size_in_morgan = C,
                            markers = R,
                            num_threads = 4,
                            time_points = t,
                            num_indiv_sampled = 20,
                            record_true_junctions = TRUE,
                            verbose = TRUE,
                            seed = 42)

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

  cat(num_j_true, obs_j, exp_j, "\n")

  testthat::expect_equal(obs_j, exp_j, tolerance = 0.2)


  N <- 10000
  R <- 10000
  t <- 20
  H_0 <- 0.5
  C <- 1

  vx <- sim_phased_unphased(pop_size = N,
                            freq_ancestor_1 = H_0,
                            total_runtime = t,
                            size_in_morgan = C,
                            markers = R,
                            num_threads = 1,
                            time_points = t,
                            num_indiv_sampled = 30,
                            seed = 42)

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