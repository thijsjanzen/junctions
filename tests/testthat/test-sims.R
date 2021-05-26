context("simulations")

test_that("finite chromosome", {
  num_markers <- 100 #chromosome size
  pop_size <- 100 # population size
  H_0 <- 0.5 # nolint
  C <- 1 # nolint
  max_time <- 500

  number_replicates <- 30
  v <- c()
  for (r in 1:number_replicates) {
    v2 <- sim_fin_chrom(pop_size, H_0, max_time, C, r, num_markers)
    v <- rbind(v, as.numeric(v2$avgJunctions) )  # nolint
  }
  v <- colMeans(v)

  predicted <- number_of_junctions(N = pop_size, R = num_markers,
                           H_0 = H_0, C = C,
                           0:max_time)

  t_sample <- sample(1:max_time, size = 10)

  for (i in t_sample) {
    testthat::expect_equal(v[i], predicted[i], tolerance = 1)
  }
})

test_that("infinite chromosome", {
  pop_size <- 50 #population size
  H_0 <- 0.5 # nolint
  C <- 1 # nolint
  max_time <- 500

  number_replicates <- 30
  v <- c()
  for (r in 1:number_replicates) {
    v2 <- sim_inf_chrom(pop_size, H_0, max_time, C, -1, r)
    v <- rbind(v, as.numeric(v2$avgJunctions) )  # nolint
  }
  v <- colMeans(v)

  predicted <- number_of_junctions(N = pop_size, R = Inf,
                                   H_0 = H_0, C = C,
                                   0:max_time)

  t_sample <- sample(1:max_time, size = 10)

  for (i in t_sample) {
    testthat::expect_equal(v[i], predicted[i], tolerance = 1)
  }

  #with 100 markers
  v <- sim_inf_chrom(pop_size, H_0, max_time, C, 100, r)  # nolint
  #with C = 20 to test Poisson > 17 code
  v <- sim_inf_chrom(10, H_0, 10, 20, -1, 42)  # nolint
  #with extremely low C, to test Poisson < 1e-6 code
  v <- sim_inf_chrom(100, H_0, 100, 1e-7, -1, 42)  # nolint
})
