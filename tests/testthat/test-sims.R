context("simulations")

test_that("finite chromosome", {
  R <- 100 #chromosome size
  N <- 100 #population size
  H_0 <- 0.5 #initial heterozygosity
  C <- 1 #number of recombinations per meiosis
  max_time <- 1000

  number_replicates <- 100
  v <- c()
  for (r in 1:number_replicates) {
    v2 <- sim_fin_chrom(N, H_0, max_time, C, r, R)
    v <- rbind(v, as.numeric(v2$avgJunctions) )  # nolint
  }
  v <- colMeans(v)

  predicted <- number_of_junctions(N = N, R = R,
                           H_0 = H_0, C = C,
                           0:max_time)

  for (i in seq_along(v)) {
    testthat::expect_equal(v[i], predicted[i], tolerance = 1)
  }
})

test_that("infinite chromosome", {

  N <- 50 #population size
  H_0 <- 0.5 #initial heterozygosity
  C <- 1 #number of recombinations per meiosis
  max_time <- 1000

  number_replicates <- 100
  v <- c()
  for (r in 1:number_replicates) {
    v2 <- sim_inf_chrom(N, H_0, max_time, C, -1, r)
    v <- rbind(v, as.numeric(v2$avgJunctions) )  # nolint
  }
  v <- colMeans(v)

  predicted <- number_of_junctions(N = N, R = Inf,
                                   H_0 = H_0, C = C,
                                   0:max_time)

  for (i in seq_along(v)) {
    testthat::expect_equal(v[i], predicted[i], tolerance = 1)
  }

  #with 100 markers
  v <- sim_inf_chrom(N, H_0, max_time, C, 100, r)  # nolint
  #with C = 20 to test Poisson > 17 code
  v <- sim_inf_chrom(10, H_0, 10, 20, -1, 42)  # nolint
  #with extremely low C, to test Poisson < 1e-6 code
  v <- sim_inf_chrom(100, H_0, 100, 1e-7, -1, 42)  # nolint
})