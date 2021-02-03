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
  ll_100 <- loglikelihood_unphased(cbind(local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   local_data$location,
                                   pop_size = 100,
                                   freq_ancestor_1 = 0.5,
                                   t = 100)

  ll_200 <- loglikelihood_unphased(cbind(local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   local_data$location,
                                   pop_size = 100,
                                   freq_ancestor_1 = 0.5,
                                   t = 200)

  testthat::expect_gte(ll_100, ll_200)

  vx <- sim_phased_unphased(pop_size = 1000,
                            freq_ancestor_1 = 0.1,
                            total_runtime = 30,
                            size_in_morgan = 1,
                            markers = 1000,
                            time_points = c(30),
                            seed = 421)

  local_data <- subset(vx, vx$individual == 0 &
                         vx$time == 30)

  ll_30 <- loglikelihood_unphased(cbind(local_data$anc_chrom_1,
                                        local_data$anc_chrom_2),
                                  local_data$location,
                                  pop_size = 1000,
                                  freq_ancestor_1 = 0.1,
                                  t = 30)

  ll_100 <- loglikelihood_unphased(cbind(local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   local_data$location,
                                   pop_size = 1000,
                                   freq_ancestor_1 = 0.1,
                                   t = 600)
  testthat::expect_gte(ll_30, ll_100)

  ll_inf <- loglikelihood_unphased(cbind(local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   local_data$location,
                                   pop_size = 1000,
                                   freq_ancestor_1 = 0.1,
                                   t = 0)
  testthat::expect_true(is.infinite(ll_inf))

  multi_ll <- loglikelihood_unphased(cbind(local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   local_data$location,
                                   pop_size = 1000,
                                   freq_ancestor_1 = 0.1,
                                   t = c(3, 3000, 300000))

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

test_that("unphased, pop size", {
  vx <- sim_phased_unphased(pop_size = 100,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 1000,
                            size_in_morgan = 1,
                            markers = 1000,
                            seed = 42)
  local_data <- subset(vx, vx$individual == 0 &
                         vx$time == 1000)
  vy <- estimate_time_unphased(cbind(local_data$anc_chrom_1,
                                     local_data$anc_chrom_2),
                               local_data$location,
                               freq_ancestor_1 = 0.5,
                               upper_lim = 2000,
                               optim_pop_size = TRUE,
                               verbose = FALSE)

  testthat::expect_equal(length(vy$par), 2)
})

test_that("cpp likelihoods", {
  vx <- sim_phased_unphased(pop_size = 100,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 1000,
                            size_in_morgan = 1,
                            markers = 1000,
                            time_points = seq(100, 900, by = 100),
                            num_threads = 1,
                            seed = 42)
  t <- 500
  local_data <- subset(vx, vx$individual == 0 &
                         vx$time == t)

  anc_matrix <- cbind(local_data$anc_chrom_1,
                      local_data$anc_chrom_2)

  for (t in c(2, 100, 1000)) {
    ll1 <- loglikelihood_unphased(local_anc_matrix = anc_matrix,
                                  locations = local_data$location,
                                  pop_size = 100,
                                  t = t,
                                  freq_ancestor_1 = 0.5)

    ll2 <- loglikelihood_unphased_cpp(local_anc_matrix = anc_matrix,
                                      locations = local_data$location,
                                      pop_size = 100,
                                      freq_ancestor_1 = 0.5,
                                      t = t,
                                      phased = FALSE)

    testthat::expect_equal(ll1, ll2)

    ll3 <- loglikelihood_phased(local_anc_matrix = anc_matrix,
                                locations = local_data$location,
                                pop_size = 100,
                                t = t,
                                freq_ancestor_1 = 0.5)

    ll4 <- loglikelihood_unphased_cpp(local_anc_matrix = anc_matrix,
                                      locations = local_data$location,
                                      pop_size = 100,
                                      freq_ancestor_1 = 0.5,
                                      t = t,
                                      phased = TRUE)

    testthat::expect_equal(ll3, ll4)
  }
})


test_that("unphased_cpp", {
  vx <- sim_phased_unphased(pop_size = 100,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 1000,
                            size_in_morgan = 1,
                            markers = 1000,
                            time_points = seq(100, 900, by = 100),
                            num_threads = 1,
                            seed = 42)

  for (t in unique(vx$time)) {
    local_data <- subset(vx, vx$individual == 0 &
                           vx$time == t)

    age1 <- estimate_time_unphased(cbind(local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   local_data$location,
                                   pop_size = 100,
                                   freq_ancestor_1 = 0.5,
                                   upper_lim = 2000,
                                   verbose = FALSE)

    age2 <- estimate_time_cpp(cbind(1,
                                             local_data$anc_chrom_1,
                                             local_data$anc_chrom_2),
                                       local_data$location,
                                       pop_size = 100,
                                       freq_ancestor_1 = 0.5,
                                       lower_lim = 2,
                                       upper_lim = 1000,
                                       verbose = FALSE,
                                       phased = FALSE)
    testthat::expect_equal(age1$minimum, age2[1], tolerance = 10)
  }
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
                            num_threads = 1,
                            time_points = 100,
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