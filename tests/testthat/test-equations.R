context("equations")

test_that("calc_K", {
  
  v1 <- calc_K(R = 10)
  expect_equal(v1, 5)
  
  v2 <- calc_K(N = 10)
  expect_equal(v2, 10)
  
  v3 <- calc_K()
  A <- is.infinite(v3)
  expect_equal(A, TRUE)
})

test_that("calculate_J and error", {
  
  maxT <- 200
  J <- calculate_J(N = 100, R = 1000, H_0 = 0.5, C = 1, maxT)
  time_estim <- calculate_time(J = J, N = 100, R = 1000, H_0 = 0.5, C = 1)
  expect_equal(time_estim, maxT)
  
  expect_silent(
    calculate_error(J = J, N = 100, R = 1000, H_0 = 0.5, C = 1, time_estim, relative = TRUE)
  )
  
  expect_silent(
    calculate_error(J = J, N = 100, R = 1000, H_0 = 0.5, C = 1, time_estim, relative = FALSE)
  )
  expect_silent(
    calculate_MAT(N = 100, R = 1000, H_0 = 0.5, C = 1)
  )
})
  
  
  


