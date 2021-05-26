context("equations")

test_that("calc_K", {
  v1 <- junctions::calc_k(R = 10)
  testthat::expect_equal(v1, 5)

  v2 <- junctions::calc_k(N = 10)
  testthat::expect_equal(v2, 10)

  v3 <- junctions::calc_k()
  a <- is.infinite(v3)
  testthat::expect_equal(a, TRUE)
})

test_that("calculate_J and error", {
  max_time <- 200
  junc <- junctions::number_of_junctions(N = 100, R = 1000,
                                      H_0 = 0.5, C = 1, max_time)

  time_estim <- junctions::estimate_time(J = junc, N = 100,
                                         R = 1000, H_0 = 0.5,
                                         C = 1)
  expect_equal(time_estim, max_time)

  expect_silent(
    junctions::time_error(N = 100, R = 1000,
                          H_0 = 0.5, C = 1, time_estim,
                          relative = TRUE)
  )

  expect_silent(
    junctions::time_error(N = 100,
                          R = 1000, H_0 = 0.5,
                          C = 1, t = time_estim,
                          relative = FALSE)
  )

  # expect NA:
  e1 <- junctions::time_error(t = c(NA),
                            N = 100,
                            R = 1000, H_0 = 0.5,
                            C = 1,
                            relative = FALSE)
  testthat::expect_true(is.na(e1))




  a1 <- junctions::time_error(t = c(5, 10, 15),
                        N = 100,
                        R = 1000, H_0 = 0.5,
                        C = 1,
                        relative = FALSE)

  a2 <- junctions::time_error(t = c(5, 10, 15),
                        N = 100,
                        R = 1000, H_0 = 0.5,
                        C = 1,
                        relative = TRUE)

  for (i in seq_along(a1)) {
    testthat::expect_lt(a2[i], a1[i])
  }



  testthat::expect_silent(
    junctions::calculate_mat(N = 100, R = 1000,
                             H_0 = 0.5, C = 1)
  )

  testthat::expect_silent(
    junctions::number_of_junctions(N = Inf, R = Inf,
                                   H_0 = 0.5, C = 1, max_time)
  )

  testthat::expect_error(
    junctions::calculate_mat(N = Inf, R = Inf,
                             H_0 = 0.5, C = 1)
  )

  testthat::expect_error(
    junctions::estimate_time(J = NA, N = 100,
                             R = 1000, H_0 = 0.5,
                             C = 1)
  )

  testthat::expect_error(
    junctions::estimate_time(N = 100, R = 1000,
                             H_0 = 0.5, C = 1)
  )

  testthat::expect_silent(
    junctions::estimate_time(J = 100, N = Inf, R = Inf,
                             H_0 = 0.5, C = 1)
  )

})
