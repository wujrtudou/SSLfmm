test_that("initialization supports equal and unequal covariance", {
  set.seed(1)
  x <- rbind(matrix(rnorm(40, -1), 20, 2), matrix(rnorm(40, 1), 20, 2))
  y <- c(rep(1, 15), rep(NA, 5), rep(2, 15), rep(NA, 5))
  a <- initialize_sslfmm(x, y, g = 2, covariance_type = "equal")
  b <- initialize_sslfmm(x, y, g = 2, covariance_type = "unequal")
  expect_equal(dim(a$sigma), c(2, 2))
  expect_equal(dim(b$sigma), c(2, 2, 2))
})

test_that("cc fit has predict method", {
  set.seed(2)
  x <- rbind(matrix(rnorm(40, -1), 20, 2), matrix(rnorm(40, 1), 20, 2))
  y <- c(rep(1, 20), rep(2, 20))
  f <- fit_sslfmm(x, y, method = "cc", covariance_type = "equal")
  pr <- predict(f, x, type = "posterior")
  expect_equal(dim(pr), c(40, 2))
  expect_equal(rowSums(pr), rep(1, 40), tolerance = 1e-7)
})

test_that("mixed known source indicator is accepted", {
  mu <- matrix(c(-1, 1), 1, 2)
  s <- simulate_mixed_missingness(80, c(.5, .5), mu, matrix(1, 1, 1),
                                  alpha = .1, mar_rate = .2, seed = 5)
  x <- as.matrix(s$data["x1"])
  f <- fit_sslfmm(x, s$data$label, g = 2, method = "mixed",
                  indicator = "known", missing_source = s$data$missing_source,
                  n_starts = 2, seed = 6)
  expect_s3_class(f, "SSLfmm")
  expect_equal(f$indicator, "known")
})

test_that("known mixed source can be supplied as logical MCAR indicator", {
  mu <- matrix(c(-1, 1), 1, 2)
  s <- simulate_mixed_missingness(80, c(.5, .5), mu, matrix(1, 1, 1),
                                  alpha = .1, mar_rate = .2, seed = 12)
  x <- as.matrix(s$data["x1"])
  f <- fit_sslfmm(x, s$data$label, g = 2, method = "mixed",
                  indicator = "known", missing_source = s$data$latent_missing,
                  n_starts = 2, seed = 13)
  expect_s3_class(f, "SSLfmm")
  expect_true(all(f$latent_missing_probability %in% c(0, 1)))
})


test_that("three-class equal-covariance MCAR direct likelihood fit is stable", {
  pi <- c(.30, .40, .30)
  mu <- matrix(c(-1.8, 0.0,
                  0.0, 1.5,
                  1.8, 0.0), nrow = 2, ncol = 3)
  sigma <- matrix(c(1.00, .25, .25, 1.00), 2, 2, byrow = TRUE)
  s <- simulate_sslfmm(
    n = 250, pi = pi, mu = mu, sigma = sigma,
    mechanism = "mcar", missing_rate = .35, seed = 123
  )
  x <- as.matrix(s$data[c("x1", "x2")])
  init <- initialize_sslfmm(x, s$data$label, g = 3,
                            covariance_type = "equal", seed = 124)
  f <- fit_sslfmm(
    x, s$data$label, g = 3, method = "mcar",
    covariance_type = "equal", init = init,
    n_starts = 3, seed = 125,
    control = list(iter.max = 1000, eval.max = 1500)
  )
  expect_s3_class(f, "SSLfmm")
  expect_true(is.finite(f$loglik))
  expect_equal(dim(f$mu), c(3, 2))
  expect_equal(dim(f$sigma), c(2, 2))
  expect_equal(rowSums(f$posterior), rep(1, nrow(x)), tolerance = 1e-7)
})

test_that("MCAR fit records optimizer and start diagnostics", {
  pi <- c(.30, .40, .30)
  mu <- matrix(c(-1.8, 0.0,
                  0.0, 1.5,
                  1.8, 0.0), nrow = 2, ncol = 3)
  sigma <- matrix(c(1.00, .25, .25, 1.00), 2, 2, byrow = TRUE)
  s <- simulate_sslfmm(
    n = 180, pi = pi, mu = mu, sigma = sigma,
    mechanism = "mcar", missing_rate = .35, seed = 321
  )
  x <- as.matrix(s$data[c("x1", "x2")])
  init <- initialize_sslfmm(x, s$data$label, g = 3,
                            covariance_type = "equal", seed = 322)
  f <- fit_sslfmm(
    x, s$data$label, g = 3, method = "mcar",
    covariance_type = "equal", init = init,
    n_starts = 3, seed = 323,
    control = list(iter.max = 1000, eval.max = 1500)
  )
  expect_true(f$optimizer %in% c("nlminb", "L-BFGS-B"))
  expect_true("optimizer" %in% names(f$start_summary))
  expect_true(is.data.frame(f$start_diagnostics))
  expect_true(all(is.finite(f$start_summary$initial_objective)))
})


test_that("g=1 entropy-dependent fitting fails informatively", {
  set.seed(123)
  x <- matrix(rnorm(30), ncol = 1)
  y <- rep(1L, 30)
  y[1:8] <- NA_integer_
  expect_error(
    fit_sslfmm(x, y, g = 1, method = "mar"),
    "posterior entropy is identically zero"
  )
  expect_error(
    fit_sslfmm(x, y, g = 1, method = "mixed"),
    "posterior entropy is identically zero"
  )
})
