test_that("simulation return contract is stable", {
  mu <- matrix(c(-1, 1), nrow = 1, ncol = 2)
  s <- simulate_mixed_missingness(
    n = 80, pi = c(.5, .5), mu = mu, sigma = matrix(1, 1, 1),
    alpha = .10, mar_rate = .25, seed = 101
  )

  expect_identical(names(s), c("data", "true_setup", "groups", "probs", "raw"))
  expect_s3_class(s$data, "data.frame")
  expect_identical(
    names(s$data)[seq_len(5 + 1L)],
    c("x1", "en", "missing", "label", "truth", "observed_missing")
  )
  expect_identical(names(s$data)[1:5], c("x1", "en", "missing", "label", "truth"))
  expect_identical(s$data$en, s$data$entropy)
  expect_identical(s$data$missing, s$data$observed_missing)
  expect_identical(s$data$missing, is.na(s$data$label))
  expect_equal(s$probs, s$data$prob_mar)

  expect_identical(names(s$true_setup)[1:3], c("pi", "mu", "sigma"))
  expect_identical(
    names(s$groups)[1:4],
    c("mar_group", "obs_group", "mcar_in_mar", "mcar_in_obs")
  )
  expect_true(all(c("observed", "mcar", "mar", "missing") %in% names(s$groups)))
  expect_setequal(c(s$groups$mar_group, s$groups$obs_group), seq_len(nrow(s$data)))
  expect_length(intersect(s$groups$mar_group, s$groups$obs_group), 0L)
  expect_setequal(c(s$groups$mcar_in_mar, s$groups$mcar_in_obs), s$groups$mcar)

  expect_true(all(c("x1", "truth", "en", "missing", "label") %in% names(s$raw)))
  expect_identical(s$raw$truth, s$data$truth)
  expect_identical(s$raw$label, s$data$label)
  expect_identical(s$raw$en, s$data$en)
})

test_that("all four simulation mechanisms share the same stable schema", {
  mu <- matrix(c(-1, 1), nrow = 1, ncol = 2)
  for (mech in c("cc", "mcar", "mar", "mixed")) {
    args <- list(
      n = 50, pi = c(.5, .5), mu = mu, sigma = matrix(1, 1, 1),
      mechanism = mech, seed = 200 + match(mech, c("cc", "mcar", "mar", "mixed"))
    )
    if (mech == "mcar") args$missing_rate <- .2
    s <- do.call(simulate_sslfmm, args)
    expect_identical(names(s), c("data", "true_setup", "groups", "probs", "raw"))
    expect_identical(names(s$data)[1:5], c("x1", "en", "missing", "label", "truth"))
    expect_identical(s$data$missing, is.na(s$data$label))
    expect_equal(s$probs, s$data$prob_mar)
  }
})

test_that("rmix has the documented simple data.frame format", {
  mu <- matrix(c(-1, 1), nrow = 1, ncol = 2)
  d <- rmix(25, c(.5, .5), mu, matrix(1, 1, 1), seed_number = 3)
  expect_s3_class(d, "data.frame")
  expect_identical(names(d), c("x1", "truth"))
  expect_equal(nrow(d), 25)
})

test_that("covariance validation gives input-level errors", {
  mu2 <- matrix(c(-1, 1, 0, 0), nrow = 2, ncol = 2)
  expect_error(
    rmix(10, c(.5, .5), mu2, matrix(1, 1, 1)),
    "dimension 2 x 2"
  )
  expect_error(
    rmix(10, c(.5, .5), mu2, matrix(c(1, 2, 0, 1), 2, 2)),
    "must be symmetric"
  )
  expect_error(
    rmix(10, c(.5, .5), mu2, matrix(c(1, 2, 2, 1), 2, 2)),
    "positive definite"
  )
  bad_arr <- array(1, c(2, 2, 3))
  expect_error(
    rmix(10, c(.5, .5), mu2, bad_arr),
    "dimension 2 x 2 x 2"
  )
})

test_that("p=1 scalar matrix and array covariance inputs are handled", {
  mu <- matrix(c(-1, 1), nrow = 1, ncol = 2)
  a <- rmix(20, c(.5, .5), mu, 1, seed_number = 1)
  b <- rmix(20, c(.5, .5), mu, matrix(1, 1, 1), seed_number = 1)
  arr <- array(1, c(1, 1, 2))
  c <- rmix(20, c(.5, .5), mu, arr, seed_number = 1)
  expect_s3_class(a, "data.frame")
  expect_s3_class(b, "data.frame")
  expect_s3_class(c, "data.frame")
})

test_that("fit output uses explicit stable aliases and class-labelled posterior", {
  mu <- matrix(c(-2, 2), nrow = 1, ncol = 2)
  s <- simulate_sslfmm(60, c(.5, .5), mu, matrix(1, 1, 1),
                        mechanism = "mcar", missing_rate = .2, seed = 301)
  x <- as.matrix(s$data["x1"])
  f <- fit_sslfmm(x, s$data$label, g = 2, method = "cc")
  expect_identical(f$en, f$entropy)
  expect_identical(f$missing, f$observed_missing)
  expect_identical(colnames(f$posterior), f$label_levels)
  expect_identical(
    names(f$start_summary)[1:7],
    c("start", "initial_objective", "objective", "convergence",
      "iterations", "optimizer", "message")
  )
})

test_that("initialization and fitting reject malformed user input clearly", {
  x <- matrix(rnorm(20), ncol = 1)
  y <- rep(c(1L, 2L), each = 10)
  expect_error(
    initialize_sslfmm(data.frame(x = x[, 1], bad = letters[1:20]), y, g = 2),
    "numeric matrix"
  )
  expect_error(
    fit_sslfmm(x, y, g = 2, method = "cc", control = 1),
    "control.*list"
  )
  bad_init <- list(pi = c(.5, NA), mu = matrix(c(-1, 1), 2, 1), sigma = matrix(1, 1, 1))
  expect_error(
    fit_sslfmm(x, y, g = 2, method = "cc", init = bad_init),
    "positive finite probabilities"
  )
})

test_that("classification posterior validation does not silently repair invalid values", {
  truth <- c(1, 2)
  pred <- c(1, 2)
  bad <- matrix(c(1.1, -0.1, .2, .8), nrow = 2, byrow = TRUE)
  colnames(bad) <- c("1", "2")
  expect_error(classification_performance(truth, pred, bad), "nonnegative")

  bad2 <- matrix(c(NA, 0, .2, .8), nrow = 2, byrow = TRUE)
  colnames(bad2) <- c("1", "2")
  expect_error(classification_performance(truth, pred, bad2), "finite numeric")
})

test_that("internal research helpers are not user API", {
  expect_false(".nll_mcar" %in% getNamespaceExports("SSLfmm"))
  expect_false(".pack_theta" %in% getNamespaceExports("SSLfmm"))
})

test_that("package version remains pinned to 0.2.2", {
  expect_identical(as.character(utils::packageVersion("SSLfmm")), "0.2.2")
})
