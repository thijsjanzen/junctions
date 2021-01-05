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

  vy <- subset(vx, vx$time == 100)

  a1 <- sum(vy$anc_chrom_1 != vy$anc_chrom_2) / length(vy$anc_chrom_1)
  expected_heterozygosity <- 2 * 0.5 * 0.5 * (1 - 1 / (2 * population_size)) ^ max_t

  testthat::expect_equal(a1, expected_heterozygosity, scale = 1, tolerance = 0.1)

  num_indiv <- length(unique(vx$individual))
  testthat::expect_equal(num_indiv, 10)
  testthat::expect_equal(length(unique(vx$time)), 2)

  found <- c()
  focal_time <- 100
  for (indiv in unique(vx$individual)) {
    local_data <- subset(vx, vx$individual == indiv &
                           vx$time == focal_time)


    estim_time <- estimate_time_unphased(cbind(local_data$anc_chrom_1,
                                               local_data$anc_chrom_2),
                                         local_data$location,
                                         pop_size = population_size,
                                         freq_ancestor_1 = 0.5,
                                         upper_lim = 200)
    found[indiv + 1] <- estim_time$minimum
  }

  testthat::expect_equal(mean(found), focal_time, tolerance = 5)

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
                                   t = c(30, 100, 300))

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

test_that("unphased_cpp", {
  vx <- sim_phased_unphased(pop_size = 100,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 1000,
                            size_in_morgan = 1,
                            markers = 1000,
                            time_points = seq(100, 900, by = 100),
                            num_threads = 1,
                            seed = 42)

  for(t in unique(vx$time)) {
    local_data <- subset(vx, vx$individual == 0 &
                           vx$time == t)

    age1 <- estimate_time_unphased(cbind(local_data$anc_chrom_1,
                                         local_data$anc_chrom_2),
                                   local_data$location,
                                   pop_size = 100,
                                   freq_ancestor_1 = 0.5,
                                   upper_lim = 2000,
                                   verbose = FALSE)

    age2 <- estimate_time_unphased_cpp(cbind(1,
                                             local_data$anc_chrom_1,
                                             local_data$anc_chrom_2),
                                       local_data$location,
                                       pop_size = 100,
                                       freq_ancestor_1 = 0.5,
                                       lower_lim = 2,
                                       upper_lim = 1000,
                                       verbose = FALSE)
    cat(t, age1$minimum, age2[1], "\n")
    testthat::expect_equal(age1$minimum, age2[1], tolerance = 10)
  }
})

test_that("unphased, threads", {

  skip_on_cran()
  population_size <- 1000

  t1 <- Sys.time()
  vx <- sim_phased_unphased(pop_size = population_size,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 1000,
                            size_in_morgan = 1,
                            markers = seq(0, 1, length.out = 100),
                            num_threads = 1,
                            seed = 43)
  t2 <- Sys.time()
  vy <- sim_phased_unphased(pop_size = population_size,
                            freq_ancestor_1 = 0.5,
                            total_runtime = 1000,
                            size_in_morgan = 1,
                            markers = seq(0, 1, length.out = 100),
                            num_threads = -1,
                            seed = 43)
  t3 <- Sys.time()

  time_one_thread = difftime(t2, t1, units = "secs")[[1]]
  time_all_threads = difftime(t3, t2, units = "secs")[[1]]
  testthat::expect_lt(time_all_threads, time_one_thread)
})
