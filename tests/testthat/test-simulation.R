test_that("p=1 covariance representations are accepted", {
  mu <- matrix(0, 1, 1)
  a <- simulate_mixed_missingness(n = 30, pi = 1, mu = mu, sigma = 1, seed = 1)
  b <- simulate_mixed_missingness(n = 30, pi = 1, mu = mu, sigma = matrix(1, 1, 1), seed = 1)
  c <- simulate_mixed_missingness(n = 30, pi = 1, mu = mu, sigma = array(1, c(1, 1, 1)), seed = 1)
  expect_true(is.list(a))
  expect_true(is.list(b))
  expect_true(is.list(c))
  expect_s3_class(a$data, "data.frame")
  expect_true(all(c("observed_missing", "latent_missing", "missing_source") %in% names(a$data)))
})

test_that("shared and component-specific covariance simulation work", {
  mu <- matrix(c(-1, 1, 0, 0), nrow = 2, ncol = 2)
  shared <- diag(2)
  unequal <- array(0, c(2, 2, 2))
  unequal[, , 1] <- diag(c(1, 2))
  unequal[, , 2] <- diag(c(2, 1))
  s1 <- simulate_sslfmm(40, c(.5, .5), mu, shared, mechanism = "mcar", missing_rate = .2, seed = 1)
  s2 <- simulate_sslfmm(40, c(.5, .5), mu, unequal, mechanism = "mar", seed = 2)
  expect_equal(s1$true_setup$covariance_type, "equal")
  expect_equal(s2$true_setup$covariance_type, "unequal")
})

test_that("one-component multivariate array covariance does not drop dimensions", {
  s <- simulate_mixed_missingness(
    n = 25, pi = 1,
    mu = matrix(0, nrow = 2, ncol = 1),
    sigma = array(diag(2), c(2, 2, 1)), seed = 11
  )
  expect_true(is.list(s))
  expect_s3_class(s$data, "data.frame")
  expect_true(all(c("x1", "x2") %in% names(s$data)))
})


test_that("g=1 constant entropy calibrates requested MAR rate", {
  target <- 0.25
  cases <- list(
    1,
    matrix(1, 1, 1),
    array(1, c(1, 1, 1))
  )
  for (S in cases) {
    s <- simulate_mixed_missingness(
      n = 60, pi = 1, mu = matrix(0, 1, 1), sigma = S,
      alpha = 0.10, mar_rate = target, xi1 = 3, seed = 99
    )
    expect_true(is.list(s))
    expect_true(is.finite(s$true_setup$xi[["xi0"]]))
    expect_equal(mean(s$probs), target, tolerance = 1e-10)
    expect_equal(unique(s$data$prob_mar), target, tolerance = 1e-10)
  }
})
